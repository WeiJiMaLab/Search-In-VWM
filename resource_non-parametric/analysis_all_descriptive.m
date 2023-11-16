%% prop corr with set size [Figure 1]
clear all; close all;

load('alldata.mat')
Nsubj = size(alldata,2);
%'params_all','nll', 'prob_corr_pred_all'
 
curr_dir = pwd;
model_fits_name = '/model_fits_cluster2/';
dirname=[curr_dir,model_fits_name];
filez=dir([curr_dir,model_fits_name]);
flst={filez.name};



Nmodels = 1;
ri_len = nan(Nsubj,Nmodels);
for si= 1:Nsubj
    for mi=1:Nmodels
        strtofind  =['model_fits_all_model_',num2str(mi),'_sbj_',num2str(si),'_run_'];
        ix=regexp(flst,strtofind);
        ix=~cellfun('isempty',ix);
        ri_len(si,mi)=sum(ix);
    end
end

npars = [5 5 4];
N_samp = 1000;
N_dsb = 100;
N_trials = 640;
delta_s_base = linspace(0.0001,pi,N_dsb)';
nbinz = 6;
nbinz_cs = 4;

condi_n  = [40 40 40 40];
set_size = [2 4 6 8];
dist     = [1 2 3 4];
%offset   = linspace(-0.2,0.2, length(dataa));

greyy =[0.7 0.7 0.7];
color_distt = [37 52 148; 65 182 196; 161 218 180; 230 220 100]'/255;
color_distt_m = (color_distt + 2*ones(3,4))/3;
msz = 3;
msz2 = 3;
fontsz = 11;
ticksz = 0.025;


params_fit = nan(Nsubj, Nmodels, max(npars));
nll_fit = nan(Nsubj, Nmodels);
prob_corr_pred_fit = nan(Nsubj, Nmodels, N_dsb,4);
prob_corr_pred = nan(Nsubj,Nmodels,N_dsb,4);
prob_corr_pred_sz = nan(Nsubj, Nmodels, 4);


prop_corr = nan(Nsubj, Nmodels, length(set_size ));
prop_corr_sd = nan(Nsubj, Nmodels, length(set_size), length(dist)); % spatial distance
prop_corr_cd = nan(Nsubj, Nmodels,length(set_size), nbinz); % color distance
prop_corr_cd_sd= nan(Nsubj, Nmodels, length(set_size), nbinz_cs, length(dist)); % color distance



model_pred = 1;


for mi=1
    cd(dirname)
    
    for si = 1:Nsubj
        
        si
        
        if model_pred
            
            nll_runs = nan(ri_len(si,mi),1);
            params_runs = nan(ri_len(si,mi),npars(mi));
            prob_corr_pred_runs = nan(ri_len(si,mi),N_dsb,4);
            
            for ri = 1:ri_len(si,mi)
                load(['model_fits_all_model_',num2str(mi),'_sbj_',num2str(si),'_run_',num2str(ri),'.mat'])
                %params_all','nll', 'prob_corr_pred_all'
                %ideally, but we still have the old setup
                nll_runs(ri) = nll;
                params_runs(ri,1:npars(mi)) = params_all;
                prob_corr_pred_runs(ri,1:N_dsb,1:4) = prob_corr_pred_all;
        
            end
            ind_min = find(nll_runs == min(nll_runs));
            
            nll_fit(si, mi) = nll_runs(ind_min(1));
            params_fit(si,mi,1:npars(mi)) = params_runs(ind_min(1),1:npars(mi));
            prob_corr_pred_fit(si,mi,1:N_dsb,1:4) = prob_corr_pred_runs(ind_min(1), 1:N_dsb,1:4);
            
            
            % size(prob_corr_pred_all) 14   100     4
            %for i = 1: 4
            prob_corr_pred = prob_corr_pred_fit;
            %end
        end
        
        
        data = alldata(si).data;
        % get the binz for each subject - could use the getBinz function, but
        % its short enough
        binz = quantile(data.col_dist, nbinz);
        binz = [0 binz max(data.col_dist)];
        binz_all(si,:) = binz;
        binz_pos_all(si,:) =  [(binz(2:end)+ binz(1:end-1))/2];
        
        
        for i = 1:length(set_size)
            N_ind           = set_size(i);
            ind_s            = data.set_size == N_ind;
            prop_corr(si,mi,i) = sum(data.response(ind_s))/(sum(ind_s));
            delta_s_col = data.col_dist';
            
            if model_pred
                prop_corr_pred_sz(si,mi,i) = mean(prob_corr_pred(si,mi,1:N_dsb,i),3);
            end
            
            % group by color dist
            for j = 1:nbinz+1
                ind_cd = data.col_dist > binz(j) & data.col_dist < binz(j+1) ;
                ind_s_cd = ind_s & ind_cd;
                prop_corr_cd(si,mi,i,j) = ( sum(data.response(ind_s_cd == 1))/sum(ind_s_cd));
                
                if model_pred
                    ind_cd_d = delta_s_base > binz(j) & delta_s_base < binz(j+1) ;
                    prop_corr_pred_cd(si,mi,i,j) = mean(prob_corr_pred(si,mi,ind_cd_d,i),3);
                    %prop_corr_cd_pred(k,i,j) = mean(prob_corr_predD(k,i,ind_cd),3);
                end
                
            end
            
        end
    end
    cd ..
    
    % finding the delta_s values
    
    allJs = squeeze(params_fit(:,1,1:4));
    allJs= exp(allJs);
    allJs_lin = reshape(allJs, 4*Nsubj,1)
    
    % proportion correct with set size plot
    figure(1)
    set(gcf, 'Position', [100 100 700 390])
    marginsa=[0.08 0.05 0.10 0.03]; %left right bottom top
    guttera=[0.08 0.18];
    
    tight_subplot(2,3,1,mi, guttera, marginsa)
    
    if model_pred
      h_p=fill([set_size set_size(end:-1:1)],   [mean(squeeze(prop_corr_pred_sz(:,mi,:)),1) - std(squeeze(prop_corr_pred_sz(:,mi,:)),1)/sqrt(Nsubj)...
         fliplr(mean(squeeze(prop_corr_pred_sz(:,mi,:)),1) + std(squeeze(prop_corr_pred_sz(:,mi,:)),1)/sqrt(Nsubj))],greyy, 'EdgeColor', 'None'); hold on;
        
    end
    
    for i = 1: length(set_size)
        plot(set_size(i), prop_corr(:, mi,i), 'o', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', msz); hold on;
    end
    leg= errorbar(set_size, mean(squeeze(prop_corr(:,mi,:)),1), std(squeeze(prop_corr(:,mi,:)),1)/sqrt(Nsubj), '-o','MarkerSize', msz2, 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k','Color','k', 'LineWidth',2); hold on;
    
    
    box off
    xlim([set_size(1)-0.5 set_size(end)+0.5])
    ylim([0.5 1])
    xlabel('Set size', 'FontName','Helvetica', 'FontSize', fontsz)
    ylabel('Proportion correct', 'FontName','Helvetica', 'FontSize', fontsz)
    %title('Proportion Correct in relation to set size', 'FontName', 'Helvetica', 'FontSize', 14)
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [ticksz ticksz])
    set(gca, 'xTick', set_size)
    set(gca, 'xTicklabels', set_size)
    
    % proportion correct with color distance plot
    figure(1)
    tight_subplot(2,3,2,mi, guttera, marginsa)
    binz_pos = mean(binz_pos_all,1);
    
    for j = 1: length(set_size)
        
        if model_pred
           h_p=fill([binz_pos binz_pos(end:-1:1)],   [mean(squeeze(prop_corr_pred_cd(:,mi,j,:)),1)- std(squeeze(prop_corr_pred_cd(:,mi,j,:)),1)/sqrt(Nsubj)...
               fliplr(mean(squeeze(prop_corr_pred_cd(:,mi,j,:)),1)+ std(squeeze(prop_corr_pred_cd(:,mi,j,:)),1)/sqrt(Nsubj))],color_distt_m(:,j)', 'EdgeColor', 'None'); hold on;
        end
        
        errorbar(binz_pos, squeeze(mean(prop_corr_cd(:,mi,j,:),1)), squeeze(std(prop_corr_cd(:,mi,j,:),1))/sqrt(Nsubj), 'o-','Color',color_distt(:,j), 'MarkerSize', 2, 'LineWidth', 2); hold on;
        % note: squeeze(binz_pos_all(si,:)) used to be binz_pos. Okay edit?
        % and also added another squeeze to the prop_corr_cd
        
        
    end
    %set(gca, 'xtick', binz_pos(1:2:7))
    set(gca, 'xtick', [0 0.5 1 1.5 2 2.5 3])
    set(gca, 'xticklabels', {'0', '', '1', '' '2', '' '3'})
    %set(gca, 'xticklabels',num2str(binz_pos(1), '%0.2f'), num2str(binz_pos(3), '%0.2f'), num2str(binz_pos(5), '%0.2f'), num2str(binz_pos(7), '%0.2f') )
    %set(gca, 'FontSize', 14)
    %title('Proportion Correct in relation to color distance', 'FontName', 'Helvetica', 'FontSize', 14)
    xlabel('Color distance (radians)', 'FontName','Helvetica', 'FontSize', fontsz)
    ylabel('Proportion correct', 'FontName','Helvetica', 'FontSize', fontsz)
    %legend('N = 2', 'N = 4', 'N = 6', 'N = 8') %% NOTE: LEGEND NOT MATCHING UP
    %legend('{\it{N}} = 2', '{\it{N}} = 4', '{\it{N}} = 6', '{\it{N}} = 8')
    %legend boxoff
    set(gca, 'tickdir', 'out')
    set(gca,'ticklength', [ticksz ticksz])
    box off
    %set(gca, 'xTick', binz_pos)
    %set(gca, 'xTicklabels', binz_pos)
    xlim([0 pi])
    ylim([0.5 1])
    
    hold on;
end
psname = ['summary_statistic_UPDATED', '.pdf'];
%print_pdf(psname)


%% params fit VP model with optimal decision rule
figure (2)
set(gcf, 'Position', [100 100 400 240])
Jbarz = squeeze(params_fit(:,1,1:4))
plot(set_size, mean(allJs,1), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
errorbar(set_size, mean(allJs,1), std(allJs)/sqrt(Nsubj)); hold on;
xlim([1.5 8.5])
xlabel('set size', 'FontName', 'Helvetica','FontSize', 1.2*fontsz)
ylabel('Mean precision J_{bar}','FontName', 'Helvetica','FontSize', 1.2*fontsz)
set(gca, 'tickdir', 'out')
set(gca,'ticklength', [ticksz ticksz])
box off

%%
psname = 'precisions_set_size_for_VPO_model_UPDATED.pdf'
%print_pdf(psname)
%% 2 way RM ANOVA
%{

for i = 1 : 4 % set size
    for j = 1 : 7 % col dist bins
        %for k = 1 : 4 % sp dist
        Y((4* (i-1)+ j-1)*Nsubj +[1:1:Nsubj]) = squeeze(prop_corr_cd_sd(:,i,j,k));
        S((4* (i-1)+j-1)*Nsubj +[1:1:Nsubj]) = 1:1: Nsubj;
        F1((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = i*ones(1,Nsubj);
        F2((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = j*ones(1,Nsubj);
        %F3((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = k*ones(1,Nsubj);
        %end
    end
    
end

stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%}
%% 3 way RM ANOVA
%size(prop_corr_cd_sd)
%lillietest for normality
%{
for i = 1 : 4 % set size
    for j = 1 : 7 % col dist bins
        for k = 1 : 4 % sp dist
        Y((28* (i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = squeeze(prop_corr_cd_sd(:,i,j,k));
        S((28* (i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = 1:1: Nsubj;
        F1((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = i*ones(1,Nsubj);
        F2((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = j*ones(1,Nsubj);
        F3((28*(i-1)+ 4*(j-1)+k-1)*Nsubj +[1:1:Nsubj]) = k*ones(1,Nsubj);
        
        end
    end
end

%Y(Y<0.01) = 0.01;
%Y = log(Y);

X = [ Y; F1; F2; F3; S]';
%[RMAOV33] = RMAOV33(X,alpha)
[RMAOV33r] = RMAOV33(X)

%}

%%
psname ='_model_M5F.pdf'
%print_pdf(psname)

%%
psname ='data_real_sum_stats_v2F.pdf'
%print_pdf(psname)

%%
psname = 'prop_corr_with_col_dist_M0F.pdf'
%print_pdf(psname)




%% repeated measures 2 way anova 
%{

FACTNAMES = {'Set size', 'Dist'}



Y = nan(1, Nsubj* length(set_size)*length(dist));
S = nan(1, Nsubj* length(set_size)*length(dist));
F1 = nan(1, Nsubj* length(set_size)*length(dist));
F2 = nan(1, Nsubj* length(set_size)*length(dist));

for i = 1 : 4
    for j = 1:4
        Y((4* (i-1)+ j-1)*Nsubj +[1:1:Nsubj]) = squeeze(prop_corr_sd(:,i,j));
        S((4* (i-1)+ j-1)*Nsubj +[1:1:Nsubj]) = 1:1: Nsubj;
        F1((4* (i-1)+ j-1)*Nsubj +[1:1:Nsubj]) = i*ones(1,Nsubj);
        F2((4* (i-1)+ j-1)*Nsubj +[1:1:Nsubj]) = j*ones(1,Nsubj);
    end
    
end

stats = rm_anova2(Y,S,F1,F2,FACTNAMES)


%% repeated measures 3 way anova
% use RMAOV33.m

    %% AIC and BIC comparison:
    
    AIC = NaN(3,14);
    AICc = NaN(3,14);
    BIC =  NaN(3,14);
    logLik = NaN(3,14);
    for  mi = 1:3

        AIC(mi,:) = (2*npars(mi)) + 2*squeeze(nll_fit(:,mi)) ;
        AICc(mi,:) = AIC(mi,:) + (2*npars(mi) * (npars(mi)+1))/(N_trials-npars(mi)-1);
        
        BIC(mi,:) = (log(N_trials) * npars(mi) ) + 2*squeeze(nll_fit(:,mi)) ; 
        
    end
    
    
    delta_AICc = AICc - AICc(1,:);
    delta_AICc = delta_AICc(2:3,:);
    
    delta_BIC = BIC - BIC(1,:);
    delta_BIC = delta_BIC(2:3,:);
    %%
    close all;
    
    figure 
    set(gcf, 'Position', [100 100 400 300])
    
    marginsa = [0.10 0.07 0.16 0.13]; %left right bottom top
    guttera = [0.14 0.18]; % btwn cols, btwtn rows
    
    eb_w = 0.15; %errorbar width
    eb_t = 0.4;  %errorbar line thickness
    fontsz = 10;
    msz = 3; 
    
    for i = 1:2 
        
        tight_subplot(1,2,1,1, guttera, marginsa)
        plot(i*ones(1,Nsubj), delta_AICc(i,:), 'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', msz); hold on;
        bb1(i) = bar(i, mean(delta_AICc(i,:)), 'FaceAlpha', .3); hold on;
        set(bb1(i), 'FaceColor', [0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5], 'BarWidth', 0.3)
        ee1(i) = errorbar(i, mean(delta_AICc(i,:)), std(delta_AICc(i,:))/sqrt(Nsubj)); hold on;
        set(ee1(i), 'Color', 'k', 'Linewidth', 1); hold on;
        errorbar(ee1(i),eb_w,eb_t); hold on;
        ylim([-1 10])
        xlim([0.5 2.5])
        set(gca, 'xtick',[1 2])
        set(gca, 'XTickLabel',{'Model_2_-_1','Model_3_-_1'});
        box off
        set(gca,'tickdir', 'out')
        set(gca, 'ticklength', [0.04 0.04])
        text(1.2,10.5, '\Delta{AICc}', 'FontName', 'Helvetica','FontSize', 1.4*fontsz)
        
        
        tight_subplot(1,2,1,2, guttera, marginsa)
        plot(i*ones(1,Nsubj), delta_BIC(i,:), 'o', 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k', 'MarkerSize', msz); hold on;
        bb2(i) = bar(i, mean(delta_BIC(i,:)), 'FaceAlpha', 0.3); hold on;
        set(bb2(i), 'FaceColor', [0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5], 'BarWidth', 0.3)
        ee2(i) = errorbar(i, mean( delta_BIC(i,:) ), std(delta_BIC(i,:))/sqrt(Nsubj)); hold on;
        set(ee2(i), 'Color', 'k', 'Linewidth', 1); hold on;
        errorbar(ee2(i),eb_w,eb_t); hold on;
        ylim([-1 10])
        xlim([0.5 2.5])
        set(gca, 'xtick',[1 2])
        set(gca, 'XTickLabel',{'Model_2_-_1','Model_3_-_1'});
        box off
        set(gca,'tickdir', 'out')
        set(gca, 'ticklength', [0.04 0.04])
        text(1.2,10.5, '\Delta{BIC}', 'FontName', 'Helvetica','FontSize', 1.4*fontsz)
         
    end
    
    
    
    %%
    
    psname ='AIC_BIC_M1-3.pdf'
   % print_pdf(psname)
    %}


