function [V,ge,gi,fe,fi] = t5_delta(params,spfr_data)
%               t,  tre,tde,ae,me,sie,tri,tdi,ai,mi,sigi,a2,m
params = params.*[10, 10, 1, 1, 10, 10, 1, 1, 1, .1];

p.Vr = 0;
p.Ve = 65;
p.Vi = -10;

%stim parameters
width = spfr_data.width-1;
p.width = spfr_data.width;
p.stimDur = spfr_data.stimDur;
pos_vect = spfr_data.pos_vect;
stim_time = fix(spfr_data.time(spfr_data.stimIdx));

M = 1:(length(spfr_data.pos_vect));
%x = sort(spfr_data.pos_vect,'ascend');
x = min(pos_vect):max(pos_vect);

%if moving bar, we need to tag stims to end of pos_vect
if ~spfr_data.cat_flag
    d = mean(diff(pos_vect)) > 0; %if diff >1, then it's PD, if it's less than 1 it's ND
    if d && width > 0
        pos_vect = [pos_vect, pos_vect(end)+(1:(width))]; %append stim positions to PD side of pos_vect
    elseif ~d && width > 0
        pos_vect = [pos_vect(1)+((width):-1:1), pos_vect];

    end
else %if spfr, we dont want a window
    x =  [x(1)+(-width:-1), x];
end
    
Ie = zeros(length(stim_time),length(M));
Ii = zeros(length(stim_time),length(M));

%pos_vect defines leadingedge of stimulated position. x defines valid
%positions in window
for i = 1:length(pos_vect)
    %pos = pos_vect(i) - min(pos_vect) + 1;
    idx = x <= pos_vect(i) & x >= pos_vect(i)-width;

     Ie(i,idx) = 1;
     Ii(i,idx) = 1;

end
p.Ie = Ie;
p.Ii = Ii;

%excitatory
p.Tre = params(1);
p.Tde = params(2);
mue = params(3);
sige = params(4);

%inhibitory
p.Tri = params(5); 
p.Tdi = params(6);   
Ai = params(7); 
mui = params(8);   
sigi = params(9);  

%time
p.b = spfr_data.stimDur*params(10); % b = m*effDur

if ~spfr_data.cat_flag
    p.b = p.b*spfr_data.width; %if moving, then also account for width 
end


p.a2 = -1; % fix a2
        
% Mode
p.ae = exp(-(x-mue).^2 / (2*sige^2) );
p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );

p.stim_time = stim_time;
y0 = zeros((4*size(Ie,2)),1);
t_ind = find([spfr_data.stimIdx;1]);
[~,s_ind] = min(abs(spfr_data.time - (spfr_data.time(spfr_data.stimIdx) + p.stimDur)'));
s_ind = s_ind';

p.stim_time = stim_time;
if spfr_data.cat_flag
    %solve the diffeq once
    %solver for during flash
    p.stim_num = 1;
    sol = ode45(@(t,y) vm_wrap_sb(t,y,p),[spfr_data.time(1) spfr_data.time(s_ind(1))+10],y0);
    tspan = spfr_data.time(t_ind(1):s_ind(1)); % solve during bar flash
    tmp{1} = deval(sol,tspan);
    %solve after the flash
    y1 = tmp{1}(:,end); %set initial conditions for second half
    y1(1:4:length(y0),end) = p.b/p.Tre; % set delta function, which jumps to amplitude instantly
    sol = ode45(@(t,y) vm_wrap_sb(t,y,p),[spfr_data.time(s_ind(1)) spfr_data.time(t_ind(2))+10],y1);
    tspan = spfr_data.time((s_ind(1)+1):(t_ind(2)-1)); % solve after flash
    tmp{2} = deval(sol,tspan);
    ode_res = [tmp{:}];    
    
    %populate a matrix with this dynamic for each stim
    he = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
    hi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    
    he_tmp = ode_res(1:4:length(y0),:);
    fe_tmp = ode_res(2:4:length(y0),:);
    hi_tmp = ode_res(3:4:length(y0),:);
    fi_tmp = ode_res(4:4:length(y0),:);
    
    for i = 1:size(Ie,1)
        he(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(he_tmp,2)-1)) = he_tmp(logical(p.Ie(1,:)),:);
        fe(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(fe_tmp,2)-1)) = fe_tmp(logical(p.Ie(1,:)),:);
        hi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(hi_tmp,2)-1)) = hi_tmp(logical(p.Ii(1,:)),:);
        fi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fi_tmp,2)-1)) = fi_tmp(logical(p.Ii(1,:)),:);
    end
    
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe;
    gi  = p.ai*fi;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));

else
    p.Ie = [p.Ie;zeros(1,size(p.Ie,2))];
    t_ind(end) = t_ind(end)-1;
    for i = 1:length(t_ind)-1
        sol = ode45(@(t,y) vm_wrap_mov(t,y,p),[spfr_data.time(t_ind(i)) spfr_data.time(t_ind(i+1))+10],y0);
        tspan = spfr_data.time(t_ind(i):t_ind(i+1)-1); % solve during bar flash
        tmp{i} = deval(sol,tspan); %solve each bit and store laterally
        gate = (p.Ie(i,:) - p.Ie(i+1,:)) > 0;
        y0 = tmp{i}(:,end); %next bit picks up where this left off
        y0(1:4:length(y0),end) = (p.b/p.Tre).*gate' + tmp{i}(1:4:length(y0),end).*~gate'; % set delta function, which jumps to amplitude instantly only at section that relaxes at t_ind
    end

    ode_res = [tmp{:}];
    ode_res = [ode_res,ode_res(:,end)];         %THIS IS BECAUSE WE SHAVED OFF LAST TIME POINT IN FOR LOOP STITCHING, FIND BETTER WAY TO DO THIS
    he = ode_res(1:4:length(y0),:);
    fe = ode_res(2:4:length(y0),:);
    hi = ode_res(3:4:length(y0),:);
    fi = ode_res(4:4:length(y0),:);

    ge  = sum(fe'.*p.ae,2);
    gi  = sum(fi'.*p.ai,2); 
    
    V = (p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi);
        
       
end

end

%ODE45
function dydt = vm_wrap_mov(t,y,p)
    p.stimi = 0;
    p.stim_numi = 1;
    id1 = t > p.stim_time ;
    id2 = t < p.stim_time + p.stimDur;
    tmp = find(id1.*id2);
    if any(tmp)
        p.stimi = 1;
        p.stim_numi = tmp(1);
    end

    dydt = vm_solve(y,p);
end

function dydt = vm_wrap_sb(t,y,p)
    p.stimi = 0;
    p.stim_numi = p.stim_num;
    if t < p.stim_time(p.stim_num) + p.stimDur
        p.stimi = 1;
    end

    dydt = vm_solve(y,p);
end

function dydt = vm_solve(y,p)

he = y(1:4:end);
fe = y(2:4:end);
hi = y(3:4:end);
fi = y(4:4:end);


dhe = (-he')   ./ p.Tre;
dfe = (-fe' + he')                           ./ p.Tde;
dhi = (-hi' + p.stimi*p.Ii(p.stim_numi,:))   ./ p.Tri;
dfi = (-fi' + hi')                           ./ p.Tdi;

tmp = [dhe;dfe;dhi;dfi];
dydt = reshape(tmp,[],1);

end