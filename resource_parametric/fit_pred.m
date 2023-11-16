function [params_fit, nll, prob_corr_pred, loglik_all] = fit_pred(mi, setsz, response, delta_s_col, N_samp, N_dsb)

N_set    = [2 4 6 8];
% fitting the params
switch mi
    
    case 1 % VPF model - Jbar, alpha, tau
        PLB  = [log([0.01]) log([0.01]) log([0.01]) ]; % lower bound
        PUB  = [log([600]) log([5]) log([600]) ];
        LB = [-Inf -Inf -Inf];
        UB = [Inf Inf Inf];
        
    case 2  % epf model - Jbar, alpha
        %PLB  = [log([0.01]) log([0.01])]; % lower bound
        PLB  = [log([0.01]) log([0.01])]; % lower bound
        PUB  = [log([600]) log([5])];
        %LB = [-Inf -Inf];
        %UB = [Inf Inf];
        LB = PLB;
        UB = PUB;
        
    case 3 % VP model - Jbar, alpha, tau
        PLB  = [log([0.01]) log([0.01]) log([0.01]) ]; % lower bound
        PUB  = [log([600]) log([5]) log([600]) ];
        LB = [-Inf -Inf -Inf];
        UB = [Inf Inf Inf];
        
    case 4  % ep model - Jbar, alpha
        PLB  = [log([0.01]) log([0.01])]; % lower bound
        PUB  = [log([600]) log([5])];
        LB = [-Inf -Inf];
        UB = [Inf Inf];
        
        
end

nvars      = length (PLB);
start_pars = (PUB-PLB).*rand(1,nvars) + PLB   % Initial point

OPTIONS = bads('defaults');             % Default options
OPTIONS.Ninit = 10;                      % Only 2 points for initial mesh
OPTIONS.UncertaintyHandling = true;        % Activate noise handling
%OPTIONS.NoiseSize = 1.5;                  % Estimated noise magnitude
OPTIONS.SpecifyTargetNoise = true;  % We are also specifying the noise
OPTIONS.NoiseFinalSamples   = 100;

neg_loglik_Ki = nan(max(N_set)+1, 1);
params_fit_Ki = nan(max(N_set)+1, nvars + 1);

if ismember(mi, [1 2]) % vpf, epf model
    for Ki = 1:1:max(N_set)+1
        K = Ki-1;
        K = 8;
        % make sure to change this K value!!!!
        %fun        = @(pars) -sum(LL_EVPF(mi,setsz,delta_s_col,response, N_samp,[pars K]));
        fun        = @(pars) -sum(LL_EVPF(mi,setsz,delta_s_col,response, N_samp,[pars K]));
        [params_fit_e, nll_e, ~, ~] = bads(fun,start_pars,LB,UB, PLB, PUB, OPTIONS);
        
        neg_loglik_Ki(Ki) = sum(nll_e);
        params_fit_Ki(Ki,:) = [params_fit_e K];
        
    end % end of Ki loop
    
    ind_min = find(neg_loglik_Ki == min(neg_loglik_Ki));
    ind_min = ind_min(1);
    params_fit = params_fit_Ki(ind_min,:);
    nll        = neg_loglik_Ki(ind_min);
    
    [loglik_all, prob_corr_pred] = LL_EVPF(mi,setsz,delta_s_col,response, N_samp,params_fit);
else
    fun        = @(pars) -sum(LL_EVP(mi,setsz,delta_s_col,response, N_samp,pars));
    [params_fit_e, nll_e, ~, ~] = bads(fun,start_pars,LB,UB, PLB, PUB, OPTIONS);
    nll = sum(nll_e);
    params_fit = params_fit_e;
    [loglik_all, prob_corr_pred] = LL_EVP(mi,setsz,delta_s_col,response, N_samp,params_fit);
    
end


end
