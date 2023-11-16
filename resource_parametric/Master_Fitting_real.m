function Master_Fitting_real(mpr_index)
mi = 1; % 1: VPF; 2: EPF, 3: VP, 4: EP
alldata = load('alldata.mat');

Nsbj = length(alldata.alldata);
Nruns = 20;
Nsubjects = Nsbj;
par = mod(floor(([mpr_index]-0)/Nruns), Nsubjects)+1; % subject index
runi = mod(([mpr_index]-0),Nruns)+1;
rng(runi); % to ensure random initialization of starting points.

N_dsb = 300;
N_samp = 500;

i = par

model_params_all = [4, 3, 3, 2];
N_models_par = model_params_all(mi);
params_all = nan(N_models_par, 1); % at most 5 values for params.

data = alldata.alldata(i).data;
setsz       = data.set_size';
response    = data.response;
delta_s_col = data.col_dist';

[params_all(1:1:N_models_par(1)), nll,  prob_corr_pred_all, loglik_all] = fit_pred(mi,setsz, response, delta_s_col, N_samp, N_dsb);

%%
switch mi
    case 1
        savefilename=['real_fits_vpf/vpf_real_model_fits_', num2str(N_samp), 'samp_','sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
    case 2
        savefilename=['real_fits_epf/testNOW_epf_real_model_fits_', num2str(N_samp), 'samp_','sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
    case 3
        savefilename=['real_fits_vp/vp_real_model_fits_', num2str(N_samp), 'samp_','sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
    case 4
        savefilename=['real_fits_ep/ep_real_model_fits_', num2str(N_samp), 'samp_','sbj_',num2str(par),'_run_',num2str(runi),'.mat'];
end
save(savefilename, 'params_all', 'nll', 'prob_corr_pred_all', 'loglik_all', '-mat')
end



