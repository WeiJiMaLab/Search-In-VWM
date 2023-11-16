function optim_normative_sim_pow(mpr_index)

load('alldata.mat')


alldata_sim = alldata;
Nsubjects = length(alldata_sim);
Nruns = 20;%10;
mi = 1;

% find subject/run index based on input
par = mod(floor(([mpr_index]-0)/Nruns), Nsubjects)+1 % subject index
runi = mod(([mpr_index]-0),Nruns)+1;  % model fitting run index, bads needs several starting points
rng(runi);


N_samp = 600;%5;%500;


% WORKS FOR SURE.
%LB = [0.01 -Inf];
%UB = [Inf -0.01];

% bounds in LOGGED form
LB = [-Inf -Inf -Inf];
UB = [Inf Inf Inf];
PLB = [log([0.05 exp(-12)]) -3] ;
PUB = [log([600 exp(-1.0)]) 3] ;


nvars      = length (PLB);


%addpath('/scratch/as11778/bads-master/')
curr_dir = pwd;
addpath([curr_dir, '/', 'bads-master/'])
%addpath([curr_dir, '/', 'MAE/Pilot_eye/bads-master/'])



OPTIONS = bads('defaults');             % Default options
OPTIONS.Ninit = 2;                      % Only 2 points for initial mesh
OPTIONS.UncertaintyHandling = 1;        % Activate noise handling
OPTIONS.NoiseSize = 1;%2;               % Estimated noise magnitude


%cd([curr_dir, '/Norm/'])
% this can be changed..in this way.

i = par;

data = alldata_sim(i).data;
setsz       = data.set_size';
response    = data.response;
delta_s_col = data.col_dist';


% rng('shuffle');
% temp = (PUB-PLB).*rand(100,2) + PLB;
% start_pars = temp(round(100.*rand(1)), :);

start_pars = (PUB-PLB).*rand(1,nvars) + PLB;   % Initial point
fun_LL_cost = @(pars) -sum(LL_costNEW_clust_pow(mi,setsz,delta_s_col,response,N_samp,pars));

[params_fit , nll, ~, ~] = bads(fun_LL_cost,start_pars,LB,UB, PLB, PUB, OPTIONS);

[loglikk, Jbars_optim] = LL_costNEW_clust_pow(mi,setsz,delta_s_col,response,N_samp,params_fit);


savefilename=['real_fits_fminbnd_pow/', num2str(N_samp),'_samp/POW_real_norm_fits_fminbnd_',num2str(Nruns),'_starts_sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
save(savefilename, 'Jbars_optim','params_fit','nll','loglikk', 'start_pars', '-mat')

end




