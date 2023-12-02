clear all; close all;
%save(filename, 'nll_fit','nll_fit_all','params_fitted','Jbars_optim_fit', 'npars', '-mat')

%% SET BEFORE RUNNING SCRIPT (CHOOSE ONE OF THE TWO OPTIONS)
%labels = {'VP','VP Non-parametric', 'VP-RR (linear cost)', 'VP-RR (power law cost)'};
labels = {'EP','EPF', 'VP', 'VPF'};

%%

if cell2mat(labels(1)) == 'EP'
    files_to_load_1  = 'model_params_summary/nll_params_model_EP.mat';
     files_to_load_2 = 'model_params_summary/nll_params_model_EPF.mat';
     files_to_load_3 = 'model_params_summary/nll_params_model_VP.mat';
     files_to_load_4 = 'model_params_summary/nll_params_model_VPF.mat';
    comparison_index = 3; % which index/model to compare all other models with
    
else
    files_to_load_1 ='model_params_summary/nll_params_model_VP.mat';
    files_to_load_2 ='model_params_summary/nll_params_model_non-parametric-VP.mat';
    files_to_load_3 = 'model_params_summary/nll_params_model_normative_linear.mat';
    files_to_load_4 ='model_params_summary/nll_params_model_normative_power_law.mat';
    comparison_index = 1; % which index/model to compare all other models with
    
end 

% Model 1: VPF; 2: EPF, 3: VP, 4: EP
load(files_to_load_1)
Nsubj = length(nll_fit);
nll_fit_all_models(1,1:Nsubj) = nll_fit;
params_fitted_M1 = params_fitted; %params_fitted; % 
Jbars_optim_fit_M1 = Jbars_optim_fit;
npars_all(1) = npars;

load(files_to_load_2)
nll_fit_all_models(2,1:Nsubj) = nll_fit;
params_fitted_M2 = params_fitted;
Jbars_optim_fit_M2 = Jbars_optim_fit;
npars_all(2) = npars; %4 J, tau, and K

load(files_to_load_3)
nll_fit_all_models(3,1:Nsubj) = nll_fit;
params_fitted_M3 = params_fitted;
Jbars_optim_fit_M3 = Jbars_optim_fit;
npars_all(3) = npars;

load(files_to_load_4)
nll_fit_all_models(4,1:Nsubj) = nll_fit;
params_fitted_M4 = params_fitted;
Jbars_optim_fit_M4 = Jbars_optim_fit;
npars_all(4) = npars;
%%
Jbars_optim_fit_ALL(1:Nsubj,1:4, 1) = Jbars_optim_fit_M1; 
Jbars_optim_fit_ALL(1:Nsubj,1:4, 2) = Jbars_optim_fit_M2; 
Jbars_optim_fit_ALL(1:Nsubj,1:4, 3) = Jbars_optim_fit_M3;
Jbars_optim_fit_ALL(1:Nsubj,1:4, 4) = Jbars_optim_fit_M4;
%%

%%

AIC = NaN(4,Nsubj); BIC = NaN(4,Nsubj);
for mi = 1:4
    AIC(mi, 1:Nsubj) = 2*npars_all(mi)+ 2*nll_fit_all_models(mi,1:Nsubj);
    BIC(mi, 1:Nsubj) = log(640)*npars_all(mi)+ 2*nll_fit_all_models(mi,1:Nsubj);
    
end


%%

diff_models_AIC = AIC- AIC(comparison_index,:); %bsxfun(@minus,AIC(1,:), AIC ); WRONG
diff_models_BIC = BIC-BIC(comparison_index,:); %bsxfun(@minus,BIC(1,:), BIC );

sample=[];
sample2=[];
nboot = 100000;


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
set(gcf, 'Position', [100 100 830 290])
marginsa = [0.13 0.11 0.09 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.08 0.08];
%colorz = [178,132,190; 128 128 0; 175 51 5]/255;
grey_shade = 0.5;
colorz = repmat([grey_shade, grey_shade, grey_shade],4,1);
msz = 5;
msz2 = 5;
fontsz = 15;
tlen1 = 0.024;
tlen2 = 0.024;
linewi = 1.6;

tight_subplot(1,3,1,1, guttera, marginsa)
for mi = 1:4
    bar(mi, d_AIC_sum(mi), 'FaceColor', colorz(mi,:), 'EdgeColor', colorz(mi,:)); hold on;
    errorbar(mi, d_AIC_sum(mi),d_AIC_sum(mi)-bci_aic(mi,1),bci_aic(mi,2)-d_AIC_sum(mi), 'Color', 'k', 'LineWidth',linewi, 'Capsize', 0); hold on;
end

%ylim([480 580])
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', labels, 'FontName', 'Helvetica','FontSize', fontsz)
ylabel('AIC(Model) - AIC(VP)', 'FontName', 'Helvetica','FontSize', fontsz)
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylim([-50 350])
ylim([-20 200])
%xtickangle(90)

tight_subplot(1,3,1,2, guttera, marginsa)
for mi = 1:4
    bb(mi)=bar(mi, d_BIC_sum(mi),  'FaceColor', colorz(mi,:), 'EdgeColor', colorz(mi,:)); hold on;
    ee(mi) = errorbar(mi, d_BIC_sum(mi),d_BIC_sum(mi)-bci_bic(mi,1),bci_bic(mi,2)-d_AIC_sum(mi) ,'Color', 'k', 'LineWidth',linewi, 'Capsize', 0); hold on;
end
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', labels)
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('BIC(Model) - BIC(VP)', 'FontName', 'Helvetica','FontSize', fontsz)
ylim([-50 350])
ylim([-20 200])
%xtickangle(90)
  
tight_subplot(1,3,1,3, guttera, marginsa)

non_parametric_vp = load('model_params_summary/nll_params_model_EP.mat');
Jbars_optim_fit_nonparametric = exp(Jbars_optim_fit); % non-parametric j's were stored in log form

set_size = [2,4,6,8];
for j = 1: length(set_size)
    if cell2mat(labels(1)) == 'EP'
        plot(set_size(j), mean(Jbars_optim_fit_M3(:,j)), 'o-', 'MarkerFaceColor','k' , 'MarkerEdgeColor', 'k'); hold on;
        errorbar(set_size(j),mean(Jbars_optim_fit_M3(:,j)), std(squeeze(Jbars_optim_fit_M3(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', 'k', 'Linewidth',linewi, 'Capsize', 0); hold on;

        plot(set_size(j), mean(Jbars_optim_fit_nonparametric(:,j)), 'o-', 'MarkerFaceColor','r' , 'MarkerEdgeColor', 'r'); hold on;
        errorbar(set_size(j),mean(Jbars_optim_fit_nonparametric(:,j)), std(squeeze(Jbars_optim_fit_nonparametric(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', 'r', 'Linewidth',linewi, 'Capsize', 0); hold on;
    else
        if j == 1
            Jbars_optim_fit_M2 = exp(Jbars_optim_fit_M2); 
        end
        plot(set_size(j), mean(Jbars_optim_fit_M1(:,j)), 'o-', 'MarkerFaceColor','k' , 'MarkerEdgeColor', 'k'); hold on;
        errorbar(set_size(j),mean(Jbars_optim_fit_M1(:,j)), std(squeeze(Jbars_optim_fit_M1(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', 'k', 'Linewidth',linewi, 'Capsize', 0); hold on;

        plot(set_size(j), mean(Jbars_optim_fit_M2(:,j)), 'o-', 'MarkerFaceColor','r' , 'MarkerEdgeColor', 'r'); hold on;
        errorbar(set_size(j),mean(Jbars_optim_fit_M2(:,j)), std(squeeze(Jbars_optim_fit_M2(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', 'r', 'Linewidth',linewi, 'Capsize', 0); hold on;
        
        plot(set_size(j), mean(Jbars_optim_fit_M4(:,j)), 'o-', 'MarkerFaceColor','b' , 'MarkerEdgeColor', 'b'); hold on;
        errorbar(set_size(j),mean(Jbars_optim_fit_M4(:,j)), std(squeeze(Jbars_optim_fit_M4(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', 'b', 'Linewidth',linewi, 'Capsize', 0); hold on;

    end
end

if cell2mat(labels(1)) == 'EP'
    plot(set_size, mean(Jbars_optim_fit_M3,1), 'k', 'Linewidth',linewi); hold on;
    plot(set_size, mean(Jbars_optim_fit_nonparametric,1), 'r', 'Linewidth',linewi); hold on;
else
    plot(set_size, mean(Jbars_optim_fit_M1,1), 'k', 'Linewidth',linewi); hold on;
    plot(set_size, mean(Jbars_optim_fit_M2,1), 'r', 'Linewidth',linewi); hold on;
    plot(set_size, mean(Jbars_optim_fit_M4,1), 'b', 'Linewidth',linewi); hold on;
end
box off
ylabel('mean precision J','FontName', 'Helvetica','FontSize', fontsz)
xlabel('Set size', 'FontName','Helvetica', 'FontSize', fontsz)
set(gca, 'tickdir', 'out')
set(gca, 'xTick', set_size)
set(gca, 'xTicklabels', set_size)
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
xlim([set_size(1)-0.5 set_size(end)+0.5])
ylim_max = 35;
ylim([-2 ylim_max])
ylim([-2 45])
text(-0.1, 50, 'C', 'FontName', 'Helvetica', 'FontSize', 1.1*fontsz)

fontsz = 10;
delta_y = 20;
x_first = 0.05;%1.7;
y_first = 190;
%t
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
legend([bb(1), bb(2), bb(3),  bb(4)], labels)
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

