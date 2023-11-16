function [params_fit nll prob_corr_pred] = fit_pred_descriptive(mi, setsz, response, delta_s_col, N_samp, N_dsb)

rng;

N_set    = [2 4 6 8];

% fitting the params
switch mi
   
    case 1 % varibable precision: Jbar np - 4 pars and tau- optimal decision rule
        LB  = [log([0.01])*ones(1,length(N_set))  log([0.01]) ]; % lower bound
        UB  = [log([600])*ones(1,length(N_set))   log([600]) ];
        PLB = LB;
        PUB = UB;    
        
    case 2  % variable precision: Jbar np- 4 pars and tau- suboptimal decision rule
        LB  = [log([0.01])*ones(1,length(N_set))  log([0.01]) ]; % lower bound
        UB  = [log([600])*ones(1,length(N_set))   log([600]) ];
        PLB = LB;
        PUB = UB; 
        
    case 3 %  fixed precision: Jbar np- 4 pars - optimal decision rule
        LB  = [log([0.01])*ones(1,length(N_set)) ];
        UB  = [log([600])*ones(1,length(N_set)) ];
        PLB = LB;
        PUB = UB;
        
end

nvars      = length (PLB);
fun        = @(pars) -sum(LL_descriptive(mi,setsz,delta_s_col,response, N_samp,pars));
start_pars = (PUB-PLB).*rand(1,nvars) + PLB;   % Initial point

OPTIONS = bads('defaults');             % Default options
OPTIONS.Ninit = 2;                      % Only 2 points for initial mesh
OPTIONS.UncertaintyHandling = 1;        % Activate noise handling
OPTIONS.NoiseSize = 1;                  % Estimated noise magnitude

[params_fit, nll, exitflag, output] = bads(fun,start_pars,LB,UB, PLB, PUB, OPTIONS);

% define delta_s for model fitting
delta_s_base = linspace(0.0001,pi,N_dsb)';
    for N_ind    = 1:length(N_set)
    
        N = N_set(N_ind);
        prob_corr_pred(:,N_ind) =  Pred_change_descriptive(mi,N,delta_s_base,params_fit, N_samp);
    end

end
