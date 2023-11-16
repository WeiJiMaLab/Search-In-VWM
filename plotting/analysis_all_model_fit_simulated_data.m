%% prop corr with set size [Figure 1]
close all; clear all;

%%%%%% IMPORTANT TO CHECK BEFORE RUNNING %%%%%%%
alldata = load('../slot-model/alldata_sim_EPF.mat');
model_pred = 1; % if you want to plot model predictions
N_samp = 5000;
mi = 2 % 1: VPF; 2: EPF, 3: VP, 4: EP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allJs_log = [0 0 0 0]
addpath('..')
addpath('../slot-model')
curr_dir = pwd;
addpath(curr_dir);

switch mi
    case 1 % VPF
        npars = 4;
        run_length = 20;
    case 2 % EPF
        npars = 3;
        run_length = 20;
    case 3 % VP
        npars = 3;
        run_length = 20;
    case 4 % EP
        npars = 2;
        run_length = 20;
end

Nsubj = size(alldata.alldata,2);
N_set = [2, 4, 6, 8];
Nmodels = 1;
ri_len = nan(Nsubj,Nmodels);

%%

N_dsb = 100;
N_trials = 640;
delta_s_base = linspace(0.0001,pi,N_dsb)';
nbinz = 5;
nbinz_cs = 4;

condi_n  = [40 40 40 40];
set_size = [2 4 6 8];
dist     = [1 2 3 4];

greyy =[0.7 0.7 0.7];
color_distt = [37 52 148; 65 182 196; 161 218 180; 230 220 100]'/255;
color_distt_m = (color_distt + 2*ones(3,4))/3;
msz = 3;
msz2 = 3;
fontsz = 11;
ticksz = 0.025;

prob_corr_pred = nan(Nsubj,N_trials);
prob_corr_pred_sz = nan(Nsubj,  4);

prop_corr = nan(Nsubj,  length(set_size ));
prop_corr_sd = nan(Nsubj, length(set_size), length(dist)); % spatial distance
prop_corr_cd = nan(Nsubj, length(set_size), nbinz); % color distance
prop_corr_cd_sd= nan(Nsubj,  length(set_size), nbinz_cs, length(dist)); % color distance

%%

%if model_pred

nll_fit               = nan(Nsubj,1);
params_fitted         = nan(Nsubj,npars);
prob_corr_pred_fit    = nan(Nsubj,N_trials);
Jbars_optim_fit       = nan(Nsubj,length(set_size));
K_optim_fit       = nan(Nsubj,length(set_size));
cost_overall_tot_fit  = nan(Nsubj,length(set_size));
start_pars_fit        = nan(Nsubj, npars);
nll_all_runs = nan(Nsubj, run_length,1);

nll_fit_all           = nan(Nsubj,run_length);
params_fitted_all     = nan(run_length, npars);
Jbars_optim_fit_all   = nan(run_length, length(set_size));
K_optim_fit_all       = nan(Nsubj, run_length, length(set_size));

start_pars_all        = nan(run_length, npars);
nll_all_subj          = nan(Nsubj, 9);
prop_corr_pred_sz = nan(Nsubj, length(N_set));
allJs_log = nan(Nsubj, length(N_set));

for si = 1:Nsubj
%     if si == 2
%         continue
%     end
%     if si == 3
%         continue
%     end
    switch mi
        case 1
            cd ('../slot-model/sim_fits_vpf/')
        case 2
            cd ('../slot-model/sim_fits_epf/')
        case 3
            cd ('../slot-model/sim_fits_vp/')
        case 4
            cd ('../slot-model/sim_fits_ep_vp_data/')
    end
    for ri = 1:run_length
        
        switch mi
            case 1
                file = ['vpf_sim_model_fits_5000samp_sbj_',num2str(si),'_run_',num2str(ri),'.mat'];
                
                if isfile(fullfile(cd, file))
                    load(file)
                else
                    continue
                    
                end
            case 2
                file = ['epf_sim_model_fits_5000samp_sbj_',num2str(si),'_run_',num2str(ri),'.mat'];
                
                if isfile(fullfile(cd, file))
                    load(file)
                else
                    continue
                    
                end
            case 3
                file = ['vp_sim_model_fits_5000samp_sbj_',num2str(si),'_run_',num2str(ri),'.mat'];
                
                if isfile(fullfile(cd, file))
                    load(file)
                else
                    continue
                    
                end
            case 4
                file = ['ep_sim_model_fits_vp_data_5000samp_sbj_',num2str(si),'_run_',num2str(ri),'.mat'];
                
                if isfile(fullfile(cd, file))
                    load(file)
                else
                    continue
                    
                end
        end
        nll_all_runs(si,ri,:)          = nll;
        nll_fit_all(si,ri)             = -sum((loglik_all));
        params_fitted_all(ri,:)        = params_all;
    end
    
    ind_min = find(nll_fit_all(si,:) == min(nll_fit_all(si,:)));
    nll_fit(si) = nll_fit_all(si,ind_min(1));
    nll_all_subj(si,:) = squeeze(nll_all_runs(si, ind_min(1),:))';
    params_fitted(si,:) = params_fitted_all(ind_min(1), :);
    
    cd ../../plotting
    
    for nind = 1:length(N_set)
        n           = N_set(nind);
        ind_s       = alldata.alldata(si).data.set_size == n;
        delta_s_col = alldata.alldata(si).data.col_dist(ind_s)';
        
        Jbar_sz = exp(params_fitted(si,1));
        alpha = exp(params_fitted(si,2));
        Jbar = Jbar_sz*n^(-alpha);
        
        switch mi
            case 1
                tau =  exp(params_fitted(si,3));
                K = params_fitted(si,4);
                params_set = [Jbar tau K];
            case 2
                K = params_fitted(si,3);
                params_set = [Jbar K];
            case 3
                tau =  exp(params_fitted(si,3));
                params_set = [Jbar tau];
            case 4
                params_set = Jbar;
        end
        
        allJs_log(si, nind) = Jbar;
        
        switch mi
            case 1
                pc = squeeze(calc_prob_corr_EVPF(delta_s_col,mi,params_set, N_samp, n))';
            case 2
                pc = squeeze(calc_prob_corr_EVPF(delta_s_col,mi,params_set, N_samp, n))';
            case 3
                pc = squeeze(calc_prob_corr_EVP(delta_s_col,mi,params_set, N_samp))';
            case 4
                pc = squeeze(calc_prob_corr_EVP(delta_s_col,mi,params_set, N_samp))';
        end
        
        pc(pc == 1)    =  1 - 1/N_samp;
        pc(pc == 0)    =  1/N_samp;
        prob_corr_pred_fit(si, ind_s) = pc;
        
        resp = alldata.alldata(si).data.response(ind_s)';
        loglik_new(si,ind_s) = resp.* log(pc)+ ...
            (1-resp).* log(1-pc);
    end
    
    data = alldata.alldata(si).data;
    binz = quantile(data.col_dist, nbinz);
    binz = [0 binz max(data.col_dist)];
    binz_all(si,:) = binz;
    binz_pos_all(si,:) =  [(binz(2:end)+ binz(1:end-1))/2];
    
    for i = 1:length(set_size)
        N_ind           = set_size(i);
        ind_s           = data.set_size == N_ind;
        prop_corr(si,i) = sum(data.response(ind_s))/(sum(ind_s));
        delta_s_col = data.col_dist';
        
        if model_pred
            prop_corr_pred_sz(si,i) = mean(prob_corr_pred_fit(si,ind_s));
        end
        
        % group by color dist
        for j = 1:nbinz+1
            ind_cd = data.col_dist > binz(j) & data.col_dist < binz(j+1) ;
            ind_s_cd = ind_s & ind_cd;
            prop_corr_cd(si,i,j) = (sum(data.response(ind_s_cd == 1))/sum(ind_s_cd));
            
            if model_pred
                prop_corr_pred_cd(si,i,j) = mean(prob_corr_pred_fit(si,ind_s_cd),2);
            else
                prop_corr_pred_cd(si,i,j) = nan;
            end
            
        end
        
    end
end
%end

%%
allJs_log = log(allJs_log);
%%
% proportion correct with set size plot
%figure(1)
close all;
figure(1)
set(gcf, 'Position', [100 100 750 240])

fontsz = 10;
tlen1 = 0.024;
tlen2 = 0.024;
linewi = 1.75;

marginsa=[0.08 0.05 0.17 0.09]; %left right bottom top
guttera=[0.11 0.18];

tight_subplot(1,3,1,1, guttera, marginsa)

if model_pred
    h_p=fill([set_size set_size(end:-1:1)],   [mean(squeeze(prop_corr_pred_sz),1) - std(squeeze(prop_corr_pred_sz),1)/sqrt(Nsubj)...
        fliplr(mean(squeeze(prop_corr_pred_sz),1) + std(squeeze(prop_corr_pred_sz),1)/sqrt(Nsubj))],greyy, 'EdgeColor', 'None'); hold on;
    
end

for i = 1: length(set_size)
    plot(set_size(i), prop_corr(:, i), 'o', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', msz); hold on;
end
leg= errorbar(set_size, mean(squeeze(prop_corr),1), std(squeeze(prop_corr),1)/sqrt(Nsubj), '-o','MarkerSize', msz2, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k','Color','k','Linewidth',linewi, 'Capsize', 0); hold on;


box off
xlim([set_size(1)-0.5 set_size(end)+0.5])
ylim([0.5 1.02])
xlabel('Set size', 'FontName','Helvetica', 'FontSize', fontsz)
ylabel('Proportion correct', 'FontName','Helvetica', 'FontSize', fontsz)
set(gca, 'tickdir', 'out')
set(gca, 'xTick', set_size)
set(gca, 'xTicklabels', set_size)
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
text(-0.5, 1, 'A', 'FontName', 'Helvetica', 'FontSize', 1.1*fontsz)


if model_pred == 0
    for si = 1:Nsubj
        data = alldata.alldata(si).data;
        binz = quantile(data.col_dist, nbinz);
        binz = [0 binz max(data.col_dist)];
        binz_all(si,:) = binz;
        binz_pos_all(si,:) =  [(binz(2:end)+ binz(1:end-1))/2];
    end
end

% proportion correct with color distance plot
tight_subplot(1,3,1,2, guttera, marginsa)
binz_pos = mean(binz_pos_all,1);

for j = 1: length(set_size)
    
    if model_pred
        h_p=fill([binz_pos binz_pos(end:-1:1)],   [mean(squeeze(prop_corr_pred_cd(:,j,:)),1)- std(squeeze(prop_corr_pred_cd(:,j,:)),1)/sqrt(Nsubj)...
            fliplr(mean(squeeze(prop_corr_pred_cd(:,j,:)),1)+ std(squeeze(prop_corr_pred_cd(:,j,:)),1)/sqrt(Nsubj))],color_distt_m(:,j)', 'EdgeColor', 'None'); hold on;
    end
    
    errorbar(binz_pos, mean(squeeze(prop_corr_cd(:,j,:)),1), squeeze(std(prop_corr_cd(:,j,:),1))/sqrt(Nsubj)', 'o-','Color',color_distt(:,j), 'MarkerSize', 2, 'Linewidth',linewi, 'Capsize', 0); hold on;
    
end
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick', [0 0.5 1 1.5 2 2.5 3])
set(gca, 'xticklabels', {'0', '', '1', '' '2', '' '3'})
xlabel('Color distance (radians)', 'FontName','Helvetica', 'FontSize', fontsz)
ylabel('Proportion correct', 'FontName','Helvetica', 'FontSize', fontsz)
set(gca, 'FontSize', fontsz)
ylim([0.5 1.02])
set(gca, 'ticklength',[tlen1 tlen2])
text(-1.0, 1, 'B', 'FontName', 'Helvetica', 'FontSize', 1.1*fontsz)


tight_subplot(1,3,1,3, guttera, marginsa)

for j = 1: length(set_size)
    plot(set_size(j), mean(allJs_log(:,j)), 'o-', 'MarkerFaceColor',color_distt(:,j) , 'MarkerEdgeColor', color_distt(:,j)); hold on;
    errorbar(set_size(j),mean(allJs_log(:,j)), std(squeeze(allJs_log(:,j)))/sqrt(Nsubj),'o-','MarkerSize', msz2, 'Color', color_distt(:,j), 'Linewidth',linewi, 'Capsize', 0); hold on;
    
end
plot(set_size, mean(allJs_log,1), 'k', 'Linewidth',1); hold on;
box off
ylabel(' mean precision J_{bar}','FontName', 'Helvetica','FontSize', fontsz)
xlabel('Set size', 'FontName','Helvetica', 'FontSize', fontsz)
set(gca, 'tickdir', 'out')
set(gca, 'xTick', set_size)
set(gca, 'xTicklabels', set_size)
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
xlim([set_size(1)-0.5 set_size(end)+0.5])
ylim_max = max(max(allJs_log)) + 5;%40;
ylim([0 ylim_max])
text(-0.1, ylim_max +ylim_max/12, 'C', 'FontName', 'Helvetica', 'FontSize', 1.1*fontsz)
switch mi
    case 1
        sgtitle('VPF sim model')
    case 2
        sgtitle('EPF sim model')
    case 3
        sgtitle('VP sim model')
    case 4
        sgtitle('EP sim model')
end

switch mi
    case 1
        psname = 'Plots/vpf_sim_model_fits.pdf';
    case 2
        psname = 'Plots/epf_sim_model_fits.pdf';
    case 3
        psname = 'Plots/vp_sim_model_fits.pdf';
    case 4
        psname = 'Plots/ep_sim_model_fits.pdf';
end
%%
%print_pdf(psname)