function prob_corr_pred =  Pred_change_VPF(mi,setsz,delta_s,params, N_samp)
% Input  : setsz is here a number in set [2 4 6 8]; delta_s is the color
% distance in radians; params determines Jbar.
% Output : prop_cor_pred determined on a trial by trial basis.


%J1bar = exp(params);
%Nsamp = 20;%200;

N = setsz;
%Jbar = J1bar/N

switch mi
    
    case 1
        Jbar  = exp(params(N/2)); % lucky that N/2 gives indices 1,2,3,4
    case 2
        Jbar  = exp(params(N/2));
    case 3 
        Jbar  = exp(params(N/2));
   
end

if ismember(mi, [1 2])
    Jbar_pars = [Jbar exp(params(5))];
else
    Jbar_pars = Jbar;
end


prob_corr_pred = calc_prob_corr_VPF(delta_s, mi, Jbar_pars, N_samp);
prob_corr_pred(prob_corr_pred == 1) =  1 - 1/N_samp;
prob_corr_pred(prob_corr_pred == 0) =  1/N_samp;

end
