function [prob_corr] = calc_prob_corr_EVPF(delta_s,mi, Jbar_pars,Nsamp, N)
% Input : color distance (delta_s) in between 0 and pi, parameter(s)
% Output: proportion correct according to simple model


Jbar = Jbar_pars(1);
if mi == 1
    tau = Jbar_pars(2);
    K = Jbar_pars(3);
else
    K = Jbar_pars(2);
end

if (K > N)
    K = N;
end

if mi == 1 % variable precisions Jbars, models 1
    Jbars = gamrnd(Jbar/tau, tau, length(delta_s),2, Nsamp);
    kappa2 = fisher2kappa(Jbars);
else %fixed precision Jbar, model 2
    kappa2 = fisher2kappa(Jbar*ones(length(delta_s),2, Nsamp));
end

% encoding
% stimulus_one is target color, always = 0
s            =  [zeros(length(delta_s),1) delta_s];
stims        =  repmat(s, [1, 1, Nsamp]);
x            =  qrandvm(stims, kappa2);

%% if S == 1 % when BOTH target and non-target are KEPT

% decision model
dv           =   -log(besseli(0, kappa2,1)) - kappa2 + kappa2.* cos(bsxfun(@minus,x,stims(:,1,:))); %will be all 0's
[~,dd]       = max(dv,[],2);
indi_dvs_both_zeroed = find(squeeze(sum(abs(dv),2))<0.000001);
dd(indi_dvs_both_zeroed) = randi([1 2],length(indi_dvs_both_zeroed),1);

dd = squeeze(dd);
dd = 2 - dd;  % option 1 - correct, option 2 - wrong


if length(delta_s)>1
    prob_corr_one = mean(dd,2);
else
    prob_corr_one = mean(dd);
end


%% S == 2 % when BOTH target and non-target stimuli are DROPPED.
prob_corr_two =  repmat (0.5,1, length(delta_s))';

%%  S == 3 % when target is dropped
dd1      = squeeze(dv(:,2,:)) < 0;

dd1 = squeeze(dd1);

if length(delta_s)>1
    prob_corr_three = mean(dd1,2);
else
    prob_corr_three = mean(dd1);
end


%%  S == 4 % non-target is dropped
dd2      = squeeze(dv(:,1,:)) > 0;
dd2 = squeeze(dd2);

if length(delta_s)>1
    prob_corr_four = mean(dd2,2);
else
    prob_corr_four = mean(dd2);
end

%% Final Calculation

p1 = (K/N)*( (K-1)/(N-1) ); % when BOTH target and non-target are KEPT
p2 = ((N-K)/N) * ((N-K-1)/(N-1)); % when BOTH target and non-target stimuli are DROPPED.
p3 = (1-p1-p2)/2; % when target OR non-target is dropped
prob_corr = (prob_corr_one)*p1 + (prob_corr_two)*p2 + (prob_corr_three)*(p3) + (prob_corr_four)*(p3);

end


