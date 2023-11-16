function [loglik, prob_corr]=  LL_EVPF(mi,setsz,delta_s_col,response,N_samp,params)
% Output : loglikelihood of Jbar fitting to the data according to simple
% model.
% Input  : setsz = set size of data; delta_s_col = col_dist based on data; params
% determines Jbar.

N_set    = [2 4 6 8];
J1bar    = exp(params(1));
alpha = exp(params(2));

switch mi
    case 1 % VPF
        tau = exp(params(3));
        K   = params(4);
    case 2 % EPF
        K   = params(3);
end

N_trials = length(delta_s_col); % 640;

loglik = nan(N_trials, 1);

for N_ind = 1:length(N_set)
    
    N = N_set(N_ind);
    
    
    ind     = find(setsz == N);
    delta_s = delta_s_col(ind);
    resp    = response(ind)';
    
    %Jbar  = Jbar_sz(N_ind);
    Jbar  = J1bar*min(N,K)^(-alpha);
    
    % have a switch statment for here!
    if mi == 1 % vpf
        Jbar_pars = [Jbar tau K];
    else % epf
        Jbar_pars = [Jbar K];
    end
    
    prob_corr = calc_prob_corr_EVPF(delta_s,mi, Jbar_pars, N_samp, N);
    
    if (sum(prob_corr == 1) > 0)
        %flag_1 = 1
    end

    if (sum(prob_corr == 0) > 0)
       % flag_0 = 1
    end


    prob_corr(prob_corr == 1) =  1 - 1/100000000;
    prob_corr(prob_corr == 0) = 1/100000000;
    
    loglik(ind) = resp.* log(prob_corr)+ ...
        (1-resp).* log(1-prob_corr);
    
end
end
