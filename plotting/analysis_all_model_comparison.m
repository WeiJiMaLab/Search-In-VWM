clear all; close all;
%save(filename, 'nll_fit','nll_fit_all','params_fitted','Jbars_optim_fit', 'npars', '-mat')


% Model 1: VPF; 2: EPF, 3: VP, 4: EP
load('model_params_summary/nll_params_model_VPF.mat')
Nsubj = length(nll_fit);
nll_fit_all_models(1,1:Nsubj) = nll_fit;
params_fitted_M1 = params_fitted(:,1); %params_fitted; % 
Jbars_optim_fit_M1 = Jbars_optim_fit;
npars_all(1) = npars;

load('model_params_summary/nll_params_model_EPF.mat')
nll_fit_all_models(2,1:Nsubj) = nll_fit;
params_fitted_M2 = params_fitted;
Jbars_optim_fit_M2 = Jbars_optim_fit;
npars_all(2) = npars; %4 J, tau, and K

load('model_params_summary/nll_params_model_VP.mat')
nll_fit_all_models(3,1:Nsubj) = nll_fit;
params_fitted_M3 = params_fitted;
Jbars_optim_fit_M3 = Jbars_optim_fit;
npars_all(3) = npars;

load('model_params_summary/nll_params_model_EP.mat')
nll_fit_all_models(4,1:Nsubj) = nll_fit;
params_fitted_M4 = params_fitted;
Jbars_optim_fit_M4 = Jbars_optim_fit;
npars_all(4) = npars;
%%
Jbars_optim_fit_ALL(1:Nsubj,1:4, 1) = Jbars_optim_fit_M1; 
Jbars_optim_fit_ALL(1:Nsubj,1:4, 2) = Jbars_optim_fit_M2; 
Jbars_optim_fit_ALL(1:Nsubj,1:4, 3) = Jbars_optim_fit_M3;
Jbars_optim_fit_ALL(1:Nsubj,1:4, 4) = Jbars_optim_fit_M4;
Jbars_optim_fit_ALL = log(Jbars_optim_fit_ALL);
%%

%{
load('alldata.mat')
sbjind = 8;%5;%1;
data_s = alldata(sbjind).data;
mi = 1; % not used
setsz = data_s.set_size';
delta_s_col = data_s.col_dist';
response = data_s.response; %data_s.response';
N_samp = 600;
params = params_fitted_M3(sbjind,:);
for si = 1:16
    si
    [loglik_across_reps(si,1:640), Jbar_sz, cost_overall, prob_corr_tot] =  LL_costNEW_clust_pow(mi,setsz,delta_s_col,response,N_samp,params);
end


nll_across_reps = -sum(loglik_across_reps,2);
filename= ['nll_across_reps_M3_subject_index_',num2str(sbjind),'.mat'];
    save(filename, 'nll_across_reps','-mat')
    
    %}
%%

AIC = NaN(4,Nsubj); BIC = NaN(4,Nsubj);
for mi = 1:4
    AIC(mi, 1:Nsubj) = 2*npars_all(mi)+ 2*nll_fit_all_models(mi,1:Nsubj);
    BIC(mi, 1:Nsubj) = log(640)*npars_all(mi)+ 2*nll_fit_all_models(mi,1:Nsubj);
    
end


%%

%
%AIC & BIC --models 1, 2, 3
% parameters ---
%M1 --- Jbars_optim_fit_M1 and the tau
%M2 --- tau and lambda_alpha ---and Jbars_optim_fit_M2-- optimized by
%observer
%M3 --- tau, lambda alpha and pow---and Jbars_optim_fit_M3--- optimized by
%observer
%substract the descriptive model from either

diff_models_AIC = AIC- AIC(3,:); %bsxfun(@minus,AIC(1,:), AIC ); WRONG
diff_models_BIC = BIC-BIC(3,:); %bsxfun(@minus,BIC(1,:), BIC );

sample=[];
sample2=[];
nboot = 500000;%100;%100000;


    for rmi = 1:4
        
        d_AIC_sum(rmi) = sum(squeeze(diff_models_AIC(rmi,:)));
        d_BIC_sum(rmi) = sum(squeeze(diff_models_BIC(rmi,:)));
        for kk=1:nboot
            sample=randsample(diff_models_AIC(rmi,:),Nsubj,1);
            d_AIC_sums(kk,rmi) = sum(sample);
            
            sample2=randsample(diff_models_BIC(rmi,:),Nsubj,1);
            d_BIC_sums(kk, rmi) = sum(sample2);
        end

    end
%%
ci_bnd_low = 0.025;
ci_bnd_high = 0.975;

    for rmi = 1:4
        bci_aic(rmi,1:2) = [quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_high)];
        bci_bic(rmi,1:2) = [quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_high)];
    end

%%


figure(1)
set(gcf, 'Position', [100 100 540 270])
marginsa = [0.13 0.11 0.09 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.11 0.11];
%colorz = [178,132,190; 128 128 0; 175 51 5]/255;
grey_shade = 0.5;
colorz = repmat([grey_shade, grey_shade, grey_shade],4,1); 
%african violet for M1,
%FIND COLORZ FOR EACH MODEL

tight_subplot(1,2,1,1, guttera, marginsa)
for mi = 1:4
    bar(mi, d_AIC_sum(mi), 'FaceColor', colorz(mi,:), 'EdgeColor', colorz(mi,:)); hold on;
    errorbar(mi, d_AIC_sum(mi),d_AIC_sum(mi)-bci_aic(mi,1),bci_aic(mi,2)-d_AIC_sum(mi), 'Color', 'k', 'LineWidth',1.3); hold on;
end

%ylim([480 580])
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', {'VP','VP Non-parametric', 'VP Normative (linear cost)', 'VP Normative (power law cost)'})
ylabel('AIC(Model) - AIC(VP)')
ylim([-20 200])

tight_subplot(1,2,1,2, guttera, marginsa)
for mi = 1:4
    bb(mi)=bar(mi, d_BIC_sum(mi),  'FaceColor', colorz(mi,:), 'EdgeColor', colorz(mi,:)); hold on;
    ee(mi) = errorbar(mi, d_BIC_sum(mi),d_BIC_sum(mi)-bci_bic(mi,1),bci_bic(mi,2)-d_AIC_sum(mi) ,'Color', 'k', 'LineWidth',1.3); hold on;
end
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', {'VP','VP Non-parametric', 'VP Normative (linear cost)', 'VP Normative (power law cost)'})
ylabel('BIC(Model) - BIC(VP)')
ylim([-20 200])
  
%lg = legend([], {'M1: descriptive', 'M2: slots+ resources', 'M3:normative', 'M4:normative with power law'})
%ylim([480 580])
fontsz = 10;
delta_y = 20;
x_first = 0.05;%1.7;
y_first = 190;
%text(x_first,  y_first ,           'M1: VP Non-parametric (5 parameters)', 'FontSize', fontsz); hold on;
%text(x_first,  y_first-delta_y ,   'M2: VP Normative linear (2 parameters)', 'FontSize', fontsz); hold on;
%text(x_first,  y_first-2*delta_y , 'M3: VP Normative power law (3 parameters)', 'FontSize', fontsz); hold on;
%text(x_first,  y_first-3*delta_y , 'M4: EP  (2 parameters)', 'FontSize', fontsz); hold on;
%0.6314    0.7908    0.3467    0.2152
%}
%%
psname = 'model_comparison_220722_E4.pdf' %'model_comparison_220522.pdf' %'model_comparison_210815.pdf'
print_pdf(psname)

%% params

%Jbars_optim_fit_ALL = log(Jbars_optim_fit_ALL);
%Jbars_optim_fit_ALL = exp(Jbars_optim_fit_ALL);
figure(2)
%set(gcf, 'Position', [100 100 660 250])
set(gcf, 'Position', [100 100 720 250])
Nvec = [ 2 4 6 8];
guttera2 = [0.07 0.11];
marginsa2 = [0.08 0.08 0.15 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]

tight_subplot(1,10,1,[1 2 3 4], guttera2, marginsa2)
for mi = 1:4
    bb(mi)=plot( Nvec, mean( squeeze(Jbars_optim_fit_ALL(:,:,mi)),1), 'o', 'MarkerFaceColor',colorz(mi,:), 'MarkerEdgeColor', colorz(mi,:)); hold on;
    ee(mi) = errorbar(Nvec, mean( squeeze(Jbars_optim_fit_ALL(:,:,mi)),1) , std(squeeze(Jbars_optim_fit_ALL(:,:,mi)),1)./sqrt(Nsubj), 'Color',  colorz(mi,:), 'LineWidth',1.3); hold on;
end
legend([bb(1), bb(2), bb(3),  bb(4)], {'descriptive', 'slot', 'normative', 'normative with power law'})
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', 2:2:8)
xlim([1.5 8.5])
xlabel('set size N')
ylabel('log Jbar') % ylabel('log Jbar')

tight_subplot(1,10,1,[5 6], guttera2, marginsa2)
bar(1, mean(params_fitted_M1),'FaceColor', colorz(1,:)); hold on;
errorbar(1, mean(params_fitted_M1), std(params_fitted_M1)/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
bar(2, mean(params_fitted_M2(:,1)),'FaceColor', colorz(2,:)); hold on;
errorbar(2, mean(params_fitted_M2(:,1)), std(params_fitted_M2(:,1))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
bar(3, mean(params_fitted_M3(:,1)),'FaceColor', colorz(3,:)); hold on;
errorbar(3, mean(params_fitted_M3(:,1)), std(params_fitted_M3(:,1))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
bar(4, mean(params_fitted_M4(:,1)),'FaceColor', colorz(4,:)); hold on;
errorbar(4, mean(params_fitted_M4(:,1)), std(params_fitted_M4(:,1))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
box off
set(gca, 'tickdir', 'out')
xlabel('Model')
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', {'M1','M2', 'M3','M4' })
ylabel('log \tau')


tight_subplot(1,10,1,[7 8], guttera2, marginsa2)

bar(1, mean(params_fitted_M3(:,2)),'FaceColor', colorz(2,:)); hold on;
errorbar(1, mean(params_fitted_M3(:,2)), std(params_fitted_M3(:,2))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
bar(2, mean(params_fitted_M4(:,2)),'FaceColor', colorz(3,:)); hold on;
errorbar(2, mean(params_fitted_M4(:,2)), std(params_fitted_M4(:,2))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
box off
set(gca, 'tickdir', 'out')
xlabel('Model')
set(gca, 'xtick', 1:1:2)
set(gca, 'xticklabels', {'M3', 'M4' })
ylabel('log \lambda_\alpha')

tight_subplot(1,10,1,[9], guttera2, marginsa2)
bar(1, mean(params_fitted_M4(:,3)),'FaceColor', colorz(3,:)); hold on;
errorbar(1, mean(params_fitted_M4(:,3)), std(params_fitted_M4(:,3))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
box off
set(gca, 'tickdir', 'out')
xlabel('Model')
set(gca, 'xtick', 1)
set(gca, 'xticklabels', {'M3' })
ylabel('power')
ylim([0 1])

tight_subplot(1,10,1,[10], guttera2, marginsa2)
bar(1, mean(K_optim_fit(:,1)),'FaceColor', colorz(3,:)); hold on;
errorbar(1, mean(K_optim_fit(:,1)), std(K_optim_fit(:,1))/sqrt(Nsubj), 'Color',  'k', 'LineWidth',1.3); hold on;
box off
set(gca, 'tickdir', 'out')
xlabel('Model')
set(gca, 'xtick', 1)
set(gca, 'xticklabels', {'M2' })
ylabel('K')
ylim([0 8])



%%
psname ='model_params_ALL_220724.pdf'  %'model_params_ALL_220524.pdf' %'model_comparison_210815.pdf'
print_pdf(psname)

