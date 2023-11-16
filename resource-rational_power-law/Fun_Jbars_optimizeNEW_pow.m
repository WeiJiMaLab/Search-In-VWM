function y = Fun_Jbars_optimizeNEW_pow(Jbar,tau,lambda_alpha,pow,mi,N_samp,N,delta_s)
% tau and lambda are exponentiated by LL_cost
C_behavioral = 1 - mean(calc_prob_corr(delta_s,mi, [Jbar tau],N_samp));
C_neural     = lambda_alpha*N* gamma(pow+Jbar/tau)/ gamma(Jbar/tau) * (tau^pow); %lambda_alpha * N * Jbar;
y = C_behavioral + C_neural;
end