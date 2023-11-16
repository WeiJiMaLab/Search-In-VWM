clear all; close all;

mi = 3; % 1: VPF; 2: EPF, 3: VP, 4: EP
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

%j_vecs = [4 3 4 4 3 5 2 6];
%alpha_vecs = [1 0.2 0.7 1.1 -0.5 0.3 0.4 0.6];
%tau_vecs = [12 1.9 -10 -4 3.4 2 -10 3];

Nsubj = length(j_vecs);

nll_all = nan(Nsubj, 1);
for si = 1: Nsubj
    Jbar_sz  = j_vecs(si);
    alpha = alpha_vecs(si);
    tau   = tau_vecs(si);
    
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
        if mi == 3  % vp model
            Jbar_pars = [Jbar tau];
        else % ep model
            Jbar_pars = Jbar;
        end
        
        prob_corr_pred = calc_prob_corr_EVP(delta_s,mi, exp(Jbar_pars), N_samp);
        
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
        case 3
            params_set = [Jbar_sz alpha tau];
        case 4
            params_set = [Jbar_sz alpha];
    end
    [loglik_all, prob_corr] = LL_EVP(mi,setsz,delta_s_base,resp, N_samp,params_set);
    nll_all(si) = sum(loglik_all);
    
    data = struct(field1, value1, field2, value2, field3, value3,  ...
        field5, value5, field6, value6, field7, value7);
    
    alldata_sim(si).data = data;
end
%%
alldata = alldata_sim;
%save_name =  ['alldata_sim_normativeE.mat']
if mi == 3
    save_name =  ['alldata_sim_VP.mat']
elseif mi == 4
    save_name =  ['alldata_sim_EP.mat']
end
%save (save_name, 'alldata')
-nll_all
