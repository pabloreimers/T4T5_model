function t5_on_spfr_fit_wrap_classic42s(task_id)

for cell_num = [19,20]
    counter = 0;
for width = [1,2,4]
for stim_dur = [40,160] %[40,160]
try
%load data to fit to
basedir = pwd;
tic    
load(sprintf('T5_spfr_structs/spfr_ds_cell%d_dur%d_width%d_val1.mat',cell_num,stim_dur,width),'spfr_ds')
fprintf('Load Time:%d\n',fix(toc))
catch
    continue
end
counter = counter+1;

%create a structure with all of the data to fit at once
spfr_field = sprintf('spfr_ds_width%d_stimDur%d',width,stim_dur);
tmp = spfr_ds.time > spfr_ds.time(spfr_ds.stimIdx)' &...
      spfr_ds.time  < spfr_ds.time(spfr_ds.stimIdx)' + spfr_ds.stimDur; % find hyperpol time
tmp = any(tmp,2); 
spfr_ds.weight_cost = ones(size(spfr_ds.baseSub)); % initalize multiplier for each time point in cost function
spfr_ds.weight_cost(tmp) = 2; %bump the weight on hyperpol section
%spfr_ds.weight_cost(tmp) = 1/sum(tmp); %set the weight of hyper and depol sections to be equal
%spfr_ds.weight_cost(~tmp) = 1/sum(~tmp);
spfr_struct.(spfr_field) = spfr_ds;
end
end

%set parameter bounds, and initialize randomly (reproducibly) within
tic

%     Tre     Tde      Ae  mue     sige       Tri      Tdi     Ai    mui     sigi     Tre2     Tde2      Ae2  mue2     sige2 
lb = [.1,      .1,    0,   -1,    0.1,       .1       .1,      0,    -5,    0.1,    .1       .1,    0,   -5,    0.1];
ub = [70,      70,    10,   7,     10,       70,      70,     20,     7,      5,   70,      70,    10,   0,     10];


res_spfr = inf;
breakout_flag = -1;
while res_spfr > 1.1 || breakout_flag > 10
    breakout_flag = breakout_flag + 1;
rng(task_id + breakout_flag*1000)
xo = unifrnd(lb,ub);
%xo(9:10) = unifrnd([-2,0.1],[1,4]);

%set objective function (least squares residual, normalized by position and weighted so that hyperpol and depol get equal weighting) and timelimit for solver
time_lim = 60*60;
obj_fun = @(param)  (sum( spfr_struct.spfr_ds_width4_stimDur160.weight_cost.*(spfr_struct.spfr_ds_width4_stimDur160.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width4_stimDur160)).^2 )/length(spfr_struct.spfr_ds_width4_stimDur160)+...
                     sum( spfr_struct.spfr_ds_width4_stimDur40.weight_cost.*(spfr_struct.spfr_ds_width4_stimDur40.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width4_stimDur40)).^2 )/length(spfr_struct.spfr_ds_width4_stimDur40));

out_fun = @(x, optimvalues, state) (isequal(state, 'interrupt') & toc > time_lim);
%c = @(x)(sum(x(5:6)) - sum(x(1:2))); %nonlinear inequality. taue > taui written in less than 0 form taui - taue < 0.
%ceq = @(x)([]);
%nonleq = @(x)deal(c(x),ceq(x));

%find the minimum of objective function within bounds
options = optimoptions('fmincon','Display','iter','OutputFcn',out_fun);
[param,residual,exitflag,output] = fmincon(obj_fun,xo,[],[],[],[],lb,ub,[],options);
fprintf('Function Time:%d\n',fix(toc))
res_spfr = 0.2*mean( abs(spfr_struct.spfr_ds_width4_stimDur160.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width4_stimDur160)))+...
                     0.2*mean( abs(spfr_struct.spfr_ds_width2_stimDur160.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width2_stimDur160)))+...
                     0.2*mean( abs(spfr_struct.spfr_ds_width1_stimDur160.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width1_stimDur160)))+...
                     0.2*mean( abs(spfr_struct.spfr_ds_width4_stimDur40.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width4_stimDur40)))+...
                     0.2*mean( abs(spfr_struct.spfr_ds_width2_stimDur40.baseSub - t5_off_model(param,spfr_struct.spfr_ds_width2_stimDur40)));

end

%store everything in a structure
x_fit.param = param;
x_fit.residual = residual;
x_fit.exitflag = exitflag;
x_fit.output =  output;
x_fit.mov = {};
x_fit.spfr = {};

%store the moving responses to save compute in post
counter = 0;
for width = [4,2,1]
for stim_dur = [160,40] %[40,160]
try
%load data to fit to
basedir = pwd;
   
load(sprintf('T5_spfr_structs/pd_ds_cell%d_dur%d_width%d_val1.mat',cell_num,stim_dur,width),'pd_ds')
load(sprintf('T5_spfr_structs/nd_ds_cell%d_dur%d_width%d_val1.mat',cell_num,stim_dur,width),'nd_ds')
load(sprintf('T5_spfr_structs/spfr_ds_cell%d_dur%d_width%d_val1.mat',cell_num,stim_dur,width),'spfr_ds')
catch
    continue
end

if width == 1 && stim_dur==40
    continue
end
counter = counter+1;
x_fit.mov{counter,1} = t5_off_model(param,pd_ds);
x_fit.mov{counter,2} = t5_off_model(param,nd_ds);
x_fit.spfr{counter} = t5_off_model(param,spfr_ds);
end
end

%store in appropriate folder
tic
foldername = [basedir,'/xfits_T5/on_classic42s/',...
    sprintf('cell%d',cell_num)];

if ~exist(foldername, 'dir')
       mkdir(foldername)
end

save([foldername,'/x_fit',num2str(task_id),'.mat'],'x_fit')
fprintf('Save time:%d\n',fix(toc))
end
end

function [V] = t5_off_model(params,spfr_data)

params = params.*[10,10,1,1,1,10,10,1,1,1,10,10,1,1,1,];

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
Ae = params(3);
mue = params(4);
sige = params(5);

%inhibitory
p.Tri = params(6); 
p.Tdi = params(7);   
Ai = params(8); 
mui = params(9);   
sigi = params(10); 

%excitatory2
p.Tre2 = params(11);
p.Tde2 = params(12);
Ae2 = params(13);
mue2 = params(14);
sige2 = params(15);


% Model
p.ae = Ae.*exp(-(x-mue).^2 / (2*sige^2) );
%p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );
p.ai = Ai.*(abs(x-mui) < sigi);
p.ae2 = Ae2.*exp(-(x-mue2).^2 / (2*sige2^2) );

p.stim_time = stim_time;
y0 = zeros((6*size(Ie,2)),1);
t_ind = find([spfr_data.stimIdx;1]);
[~,s_ind] = min(abs(spfr_data.time - (spfr_data.time(spfr_data.stimIdx) + p.stimDur)'));
s_ind = s_ind';

p.stim_time = stim_time;

if spfr_data.cat_flag
    p.stim_num = 1;
    sol = ode45(@(t,y) vm_wrap_sb(t,y,p),[spfr_data.time(1) spfr_data.time(s_ind(1))+10],y0);
    tspan = spfr_data.time(t_ind(1):s_ind(1)); 
    tmp{1} = deval(sol,tspan);
    %solve after the flash
    y1 = tmp{1}(:,end); %set initial conditions for second half
    sol = ode45(@(t,y) vm_wrap_sb(t,y,p),[spfr_data.time(s_ind(1)) spfr_data.time(t_ind(2))+10],y1);
    tspan = spfr_data.time((s_ind(1)+1):(t_ind(2)-1)); % solve after flash
    tmp{2} = deval(sol,tspan);
    ode_res = [tmp{:}];    
    
    %populate a matrix with this dynamic for each stim
    he = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
    hi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    he2 = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fe2 = zeros(size(p.Ie,2),size(spfr_data.time,1));
    
    he_tmp = ode_res(1:6:length(y0),:);
    fe_tmp = ode_res(2:6:length(y0),:);
    hi_tmp = ode_res(3:6:length(y0),:);
    fi_tmp = ode_res(4:6:length(y0),:);
    he_tmp2 = ode_res(5:6:length(y0),:);
    fe_tmp2 = ode_res(6:6:length(y0),:);
    
    for i = 1:size(Ie,1)
        he(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(he_tmp,2)-1)) = he_tmp(logical(p.Ie(1,:)),:);
        fe(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(fe_tmp,2)-1)) = fe_tmp(logical(p.Ie(1,:)),:);
        hi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(hi_tmp,2)-1)) = hi_tmp(logical(p.Ii(1,:)),:);
        fi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fi_tmp,2)-1)) = fi_tmp(logical(p.Ii(1,:)),:);
        he2(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(he_tmp,2)-1)) = he_tmp2(logical(p.Ie(1,:)),:);
        fe2(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(fe_tmp,2)-1)) = fe_tmp2(logical(p.Ie(1,:)),:);      
    end
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));
    
else
    sol = ode45(@(t,y) vm_wrap_mov(t,y,p),[min(spfr_data.time) max(spfr_data.time)],y0);
    tspan = spfr_data.time;
    tmp = deval(sol,tspan); % V is now a 
    fe = tmp(2:6:length(y0),:);
    fi = tmp(4:6:length(y0),:);
    fe2 = tmp(6:6:length(y0),:);
    
    ge  = sum(fe'.*p.ae,2) + sum(fe2'.*p.ae2,2);
    gi  = sum(fi'.*p.ai,2); 
    
    V = (p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi);
    he = tmp(1:4:length(y0),:);
    hi = tmp(3:4:length(y0),:);
    
end

end

%ODE45
function dydt = vm_wrap_sb(t,y,p)
    p.stim = 0;
    if t < p.stim_time(p.stim_num) + p.stimDur
        p.stim = 1;
    end

    dydt = vm_solve(y,p);
end

function dydt = vm_wrap_mov(t,y,p)
    p.stim = 0;
    p.stim_num = 1;
    id1 = t > p.stim_time ;
    id2 = t < p.stim_time + p.stimDur;

    tmp = find(id1.*id2);
    if any(tmp)
        p.stim = 1;
        p.stim_num = tmp(1);
    end
    dydt = vm_solve(y,p);
end


function dydt = vm_solve(y,p)
he = y(1:6:end);
fe = y(2:6:end);
hi = y(3:6:end);
fi = y(4:6:end);
he2 = y(5:6:end);
fe2 = y(6:6:end);

dhe = (-he' + p.stim*p.Ie(p.stim_num,:)) / p.Tre;
dfe = (-fe' + he')     / p.Tde;
dhi = (-hi' + p.stim*p.Ii(p.stim_num,:)) / p.Tri;
dfi = (-fi' + hi')     / p.Tdi;
dhe2 = (-he2' + p.stim*p.Ie(p.stim_num,:)) / p.Tre2;
dfe2 = (-fe2' + he2')     / p.Tde2;


tmp = [dhe;dfe;dhi;dfi;dhe2;dfe2];
dydt = reshape(tmp,[],1);

end