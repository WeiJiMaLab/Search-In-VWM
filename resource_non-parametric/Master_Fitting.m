function master_fitting(mpr_index)
load('alldata.mat');

Nsbj = length(alldata);
decision_rule = 1; % 1 = optimal, 2 = suboptimal
curr_dir = pwd;
Nruns = 20;
Nsubjects = Nsbj;
par = mod(floor(([mpr_index]-0)/Nruns), Nsubjects)+1; % subject index
runi = mod(([mpr_index]-0),Nruns)+1;

N_set = [2 4 6 8];
Nsubj = length(alldata);
N_dsb = 100;
N_samp = 1000; 
N_models = 1;
N_models_par = 5; 

params_all = nan(N_models_par, 1); % at most 5 values for params. 
prob_corr_pred_all = nan(N_dsb, 4);


nbinz = 4; %6
psname = 'each_sbj_fit.ps';

mi = 1; % VP encoding; optimal decision rule
i = par;
data = alldata(i).data;
setsz       = data.set_size';
response    = data.response;
delta_s_col = data.col_dist'; 
    
[params_all(1:1:N_models_par(mi)), nll,  prob_corr_pred_all] = fit_pred_descriptive(mi,setsz, response, delta_s_col, N_samp, N_dsb);

%%

savefilename=['real_fits_descriptive/descriptive_model_fits_', num2str(N_samp), 'samp_','sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
save(savefilename, 'params_all', 'nll', 'prob_corr_pred_all',  '-mat')
end



