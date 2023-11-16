function [loglik, Jbar_sz, cost_overall, prob_corr_tot] =  LL_costNEW_clust(mi,setsz,delta_s_col,response,N_samp,params)
% Output : loglikelihood of Jbar fitting to the data according to simple
% model.
% Input  : setsz = [2 4 6 8]; delta_s_col = col_dist based on data; params
% determines Jbar.
N_set    = [2 4 6 8];
tau = exp(params(1));
lambda_alpha = exp(params(2));

N_trials = length(delta_s_col); % 640;
mpts = 3; % SHOULD try different starting points!!!

loglik           = nan(N_trials, 1);
cost_overall     = NaN(length(N_set), 1);
Jbar_sz          = NaN(length(N_set), 1);
Jbar_sz_mi       = NaN(mpts, 1);
cost_overall_tot = NaN(mpts, 1);
prob_corr_tot    = NaN(N_trials, 1);

% bounds in EXPONENTIATED form
%PLB  = exp(log(0.01)); % lower bound
%PUB  = exp(log(100.0));
PLB = exp(log(0.000001));
PUB = exp(log(1000.0));%exp(log(10000.0));
nvars = length(PLB);

%options = optimoptions('fmincon','MaxFunEvals',100000,'TolFun', 1e-12, 'TolX',1e-12);  %'Display', 'iter');
options = optimset('MaxFunEvals',100000.0,'TolFun', 1e-16, 'TolX',1e-16);  %'Display', 'iter');


for N_ind = 1:length(N_set)% 1%:length(N_set)
    N_ind;
    N = N_set(N_ind);
    
    ind     = find(setsz == N);
    delta_s = delta_s_col(ind);
    resp    = response(ind);
    %Npts = length(delta_s);
    
    fff = @(Jbarr) Fun_Jbars_optimizeNEW_pow(Jbarr,tau,lambda_alpha,mi,N_samp,N,delta_s);
    
    for msi = 1:mpts
        rng(msi);
        
        start_pars = (PUB-PLB).*rand(1,nvars) + PLB;   % Initial point
        %[Jbar_sz_mi(msi), cost_overall_tot(msi), ~, ~] = bads(fff,start_pars,LB,UB, PLB, PUB, OPTIONS);
        %[Jbar_sz_mi(msi),cost_overall_tot(msi)] = fmincon(fff,start_pars,[],[],[],[], LB, UB,[],options);
        %[Jbar_sz_mi(msi),cost_overall_tot(msi)] = fminsearch(fff,start_pars,options);
        [Jbar_sz_mi(msi),cost_overall_tot(msi)] = fminbnd(fff,PLB,PUB,options);
    end
    ind_min = find(cost_overall_tot== min(cost_overall_tot));
    if isempty(ind_min) % if the cost_overall_tot vector is full of NaN
        ind_min = 1;
        cost_overall(N_ind) = cost_overall_tot(ind_min(1));
        Jbar_sz(N_ind) = NaN;
    else
        cost_overall(N_ind) = cost_overall_tot(ind_min(1));
        Jbar_sz(N_ind) = Jbar_sz_mi(ind_min(1));
    end
    Jbar  = Jbar_sz(N_ind);
    
    if ismember(mi , [1 2])
        Jbar_pars = [Jbar tau];
    else
        Jbar_pars = Jbar;
    end
    
    prob_corr = calc_prob_corr(delta_s,mi, Jbar_pars, N_samp);
    
    prob_corr(prob_corr == 0) =  1/N_samp; % added back 01/28/2021
    prob_corr(prob_corr == 1) =  1 - 1/N_samp;
    
    resp = resp';
    loglik(ind) = resp.* log(prob_corr)+ ...
        (1-resp).* log(1-prob_corr);
    
    %ind_s = find(setsz == N);
    prob_corr_tot(setsz == N) = prob_corr;
end

end

