clear all; close all;
addpath('..') % add path to the previous directory to 

rng(0, "twister");

mi = 2; % 1: VPF; 2: EPF, 3: VP, 4: EP
condi_n  = [40 40 40 40];
Nset = [2 4 6 8];
sp_dist     = [1 2 3 4];

Ndsb = 640;
N_trials = 640;
Nruns = 8;
N_samp = 1000;
mpts = 3;

delta_s_base = linspace(0.0001,pi,Ndsb)';
nbinz = 6;
nbinz_cs = 4;

j_range = [3 4];
a_range = [-0.2231 0.5878];%[0.8 1.8];
t_range = [1 3];
K_range = [3 8];

Nsubj = 10;
j_vecs = (j_range(2)-j_range(1)).*rand(Nsubj,1) + j_range(1);
alpha_vecs = (a_range(2)-a_range(1)).*rand(Nsubj,1) + a_range(1);
tau_vecs = (t_range(2)-t_range(1)).*rand(Nsubj,1) + t_range(1);
k_vecs = randi(K_range, Nsubj,1);

all_params = [j_vecs alpha_vecs tau_vecs]

nll_all = nan(Nsubj, 1);
for si = 1: Nsubj
    Jbar_sz  = exp(j_vecs(si));
    alpha    = exp(alpha_vecs(si));
    tau      = exp(tau_vecs(si));
    K        = k_vecs(si);
    
    all_set_sizes = [ Nset(1)*ones(1,N_trials/4) Nset(2)*ones(1,N_trials/4)...
        Nset(3)*ones(1,N_trials/4) Nset(4)*ones(1,N_trials/4)]';
    
    
    all_sp_dist = [sp_dist(1)*ones(1,N_trials/4) sp_dist(2)*ones(1,N_trials/4)...
        sp_dist(3)*ones(1,N_trials/4) sp_dist(4)*ones(1,N_trials/4)];
    
    
    field1 = 'set_size';
    perm_ind = randperm(N_trials);
    value1 = [all_set_sizes(perm_ind)]'; %[rec.set];
    
    field7 = 'col_dist';
    value7 = delta_s_base'; %pi* rand(N_trials,1);
    
    prob_corr_pred_all = NaN(N_trials, 1);
    
    for N_ind = 1:length(Nset)
        
        N = Nset(N_ind);
        
        ind     = find(value1 == N);
        delta_s = delta_s_base(ind);
        
        Jbar  = Jbar_sz*N^(-alpha);
        switch mi
            case 1 % vpf
                Jbar_pars = [Jbar tau K];
                prob_corr_pred = calc_prob_corr_EVPF(delta_s,mi, Jbar_pars, N_samp, N);
            case 2 % epf
                Jbar_pars = [Jbar K];
                prob_corr_pred = calc_prob_corr_EVPF(delta_s,mi, Jbar_pars, N_samp, N);
            case 3  % vp model
                Jbar_pars = [Jbar tau];
                prob_corr_pred = calc_prob_corr_EVP(delta_s,mi, Jbar_pars, N_samp);
            case 4 % ep model
                Jbar_pars = Jbar;
                prob_corr_pred = calc_prob_corr_EVP(delta_s,mi, Jbar_pars, N_samp);
        end
        
        prob_corr_pred(prob_corr_pred == 1) =  1 - 1/N_samp;
        prob_corr_pred(prob_corr_pred == 0) =  1/N_samp;
        
        prob_corr_pred_all(ind) = prob_corr_pred;
    end
    
    
    field2 = 'stim_order';
    value2 = nan(8,N_trials);%[rec.ord];
    
    field3 = 'spatial_dist';
    value3 = all_sp_dist(randperm(N_trials));%[rec.dis];
    
    field5 = 'response';
    value5 =  binornd(1,prob_corr_pred_all);
    
    field6 = 'reaction_time';
    value6 = nan(1,N_trials);%[rec.rt];
    
    setsz = value1;
    resp = value5';
    switch mi
        case 1
            params_set = [log(Jbar_sz) log(alpha) log(tau) K];
            [loglik_all, prob_corr] = LL_EVPF(mi,setsz,delta_s_base,resp, N_samp,params_set);
        case 2
            params_set = [log(Jbar_sz) log(alpha) K];
            [loglik_all, prob_corr] = LL_EVPF(mi,setsz,delta_s_base,resp, N_samp,params_set);
        case 3
            params_set = log([Jbar_sz alpha tau]);
            [loglik_all, prob_corr] = LL_EVP(mi,setsz,delta_s_base,resp, N_samp,params_set);
        case 4
            params_set = log([Jbar_sz alpha]);
            [loglik_all, prob_corr] = LL_EVP(mi,setsz,delta_s_base,resp, N_samp,params_set);
    end
    
    nll_all(si) = -sum(loglik_all);
    
    data = struct(field1, value1, field2, value2, field3, value3,  ...
        field5, value5, field6, value6, field7, value7);
    
    alldata_sim(si).data = data;
end
%%
alldata = alldata_sim;
%save_name =  ['alldata_sim_normativeE.mat']
switch mi 
    case 1
        save_name =  ['alldata_sim_VPF.mat']
    case 2
        save_name =  ['alldata_sim_EPF.mat']
    case 3
        save_name =  ['alldata_sim_VP.mat']
    case 4
        save_name =  ['alldata_sim_EP.mat']
end
nll_all
%%
cd ..
%save (save_name, 'alldata')
cd 'creating sim data'/

