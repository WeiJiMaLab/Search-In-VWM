function [prob_corr] = calc_prob_corr_EVP(delta_s,mi, Jbar_pars,Nsamp)
% Input : color distance (delta_s) in between 0 and pi, parameter(s)
% Output: proportion correct according to simple model

%also take into account which decision rule is being used.

prob_corr = nan(length(delta_s),1);

if length(Jbar_pars) == 1 %fixed precision Jbar, model 3
    Jbar = Jbar_pars;
    kappa2 = fisher2kappa(Jbar*ones(length(delta_s),2, Nsamp));
    
else % variable precisions Jbars, models 1 or 2
    
    Jbars = gamrnd(Jbar_pars(1)/Jbar_pars(2), Jbar_pars(2), length(delta_s),2, Nsamp);
    kappa2 = fisher2kappa(Jbars);
    
end

% encoding
% stimulus_one is target color, always = 0
s            =  [zeros(length(delta_s),1) delta_s];
stims        =  repmat(s, [1, 1, Nsamp]);
x            =  qrandvm(stims, kappa2);

% decision model

dv           = -log(besseli(0, kappa2,1)) - kappa2 + kappa2.* cos(bsxfun(@minus,x,stims(:,1,:))); %will be all 0's
[~,dd]       = max(dv,[],2);
indi_dvs_both_zeroed = find(squeeze(sum(abs(dv),2))<0.000001);
dd(indi_dvs_both_zeroed) = randi([1 2],length(indi_dvs_both_zeroed),1);

dd = squeeze(dd);
dd = 2 - dd;  % option 1 - correct, option 2 - wrong

if length(delta_s)>1
    prob_corr = mean(dd,2);
else
    prob_corr = mean(dd);
end


end


