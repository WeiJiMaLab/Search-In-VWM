function [loglik, prob_corr] =  LL_EVP(mi,setsz,delta_s_col,response,N_samp,params)
% Output : loglikelihood of Jbar fitting to the data according to simple
% model.
% Input  : setsz = [2 4 6 8]; delta_s_col = col_dist based on data; params
% determines Jbar.
N_set    = [2 4 6 8];
%J1bar    = params(1);
switch mi
    case 3
        Jbar_sz    = exp(params(1));
        alpha = exp(params(2));
        tau = exp(params(3));
    case 4
        Jbar_sz    = exp(params(1)); 
        alpha = exp(params(2));
end

N_trials = length(delta_s_col); % 640;

loglik = nan(N_trials, 1);

for N_ind = 1:length(N_set)
    
    N = N_set(N_ind);
    
    
    ind     = find(setsz == N);
    delta_s = delta_s_col(ind);
    resp    = response(ind)';
    
    Jbar  = Jbar_sz*N^(-alpha);
    if mi == 3  % vp model
        Jbar_pars = [Jbar tau];   
    else % ep model
        Jbar_pars = Jbar;
    end
    
    prob_corr = calc_prob_corr_EVP(delta_s,mi, Jbar_pars, N_samp);
    
    prob_corr(prob_corr == 1) =  1 - 1/N_samp;
    prob_corr(prob_corr == 0) =  1/N_samp;
    
    loglik(ind) = resp.* log(prob_corr)+ ...
        (1-resp).* log(1-prob_corr);
    
end
end
