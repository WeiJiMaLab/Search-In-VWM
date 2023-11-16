function [prob_corr] = calc_prob_corr(delta_s,mi,pars,Nsamp)
% Input : color distance (delta_s) in between 0 and pi, parameter(s)
% Output: proportion correct according to simple model

prob_corr = nan(length(delta_s),1);

Jbar = pars(1);

if length(pars) == 2
    tau = pars(2);
elseif length(pars) == 3
    p_drop = pars(3);
end
% kappa2fisher takes in a Jbar and a "drop". A low drop affect the
% effeciency?

if length(pars) == 1 %fixed precision Jbar, model 3
    Jbars = Jbar*ones(length(delta_s),2, Nsamp);
else % variable precisions Jbars, models 1 or 2
    Jbars = gamrnd(Jbar/tau, tau, length(delta_s),2, Nsamp);
end

kappa2 = fisher2kappa(Jbars);

% encoding
% stimulus_one is target color, always = 0
s            =  [zeros(length(delta_s),1) delta_s];
stims        =  repmat(s, [1, 1, Nsamp]);
x            =  qrandvm(stims, kappa2);

% decision model
if ismember(mi,[1 3]) % optimal decision rule
    dv           =   -log(besseli(0, kappa2,1)) - kappa2 + kappa2.* cos(bsxfun(@minus,x,stims(:,1,:))); %will be all 0's
    [~,dd]       = max(dv,[],2);
    indi_dvs_both_zeroed = find(squeeze(sum(abs(dv),2))<0.000001);
    dd(indi_dvs_both_zeroed) = randi([1 2],length(indi_dvs_both_zeroed),1);
else % suboptimal decision rule
    dv           =  abs(angle(exp(1i*x)./exp(1i*0))); %abs(circ_dist(x,0));
    [~,dd]       = min(dv,[],2);
end

dd = squeeze(dd);
dd = 2 - dd;  % option 1 - correct, option 2 - wrong

if length(delta_s)>1
    prob_corr = mean(dd,2);
else
    prob_corr = mean(dd);
end

end
