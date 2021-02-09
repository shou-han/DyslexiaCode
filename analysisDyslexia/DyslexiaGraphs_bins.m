% establish the groups 
Methods={'HAPPE'};
allnamefile={'35'};
no_of_bins=2;
for methodindx=1
    method = Methods{methodindx};
    clearvars -except method Methods  methodindx allnamefile no_of_bins
    clc
    close all
    set(0,'DefaultFigureVisible','on')
    STFT_timer=[];
    addpath(genpath('../CSDtoolbox/'));
    addpath(genpath('../eeglab13_6_5b/'));
    addpath('functions/');
    eeglab
    % 32,20 26 has approx 5%  trials accepted due to blinking for eye tracker
    % 
    
    %4,63,65,
    if strcmp(method,'eye')
        ssj{1} = [1 2 10 16 19 21 24 27 30 36 37 41 42 44 46 51   60 62 67 68];%AM 47 58 58
        ssj{2} = [8   17 18 25 49 50 55 57 59 61 69  71 72];%RM 14 14
        ssj{3} = [3 5 7 9 11 12 13 15 22 23  28 29 31 33 34  38  40 43 45 48 52   54  53 56 35 39 66 64  73 74  75];%Dysl
    else
        ssj{1} = [1 2 10 16 19 21 24 27 30 32 36 37 41 42 44 46 51   60 62 67 68];%AM 47 58 58
        ssj{2} = [6 8   17 18 20 25 49 50 55 57 59 61 69  71 72];%RM 14 14
        ssj{3} = [3 5 7 9 11 12 13 15 22 23 26 28 29 31 33 34  38  40 43 45 48 52   54  53 56 35 39 66 64  73 74  75];%Dysl
    end
    savenameSPSS = ['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix150.csv'];
    dyscAMRM = {'AM_control','RM_control','Dyslexic'};
    
    CSD = 1;
    
    count = 0;
    Idlevel= {'Subject','group',...
        'RTLSLH','RTRSLH','RTLSRH','RTRSRH',...
        'CPPonsetLSLH','CPPonsetRSLH','CPPonsetLSRH','CPPonsetRSRH',...
        'CPPslopeLSLH','CPPslopeRSLH','CPPslopeLSRH','CPPslopeRSRH',...
        'CPPampatresponseLSLH','CPPampatresponseRSLH','CPPampatresponseLSRH','CPPampatresponseRSRH',...
        'CPPPeakTimeLSLH','CPPPeakTimeRSLH','CPPPeakTimeLSRH','CPPPeakTimeRSRH',...
        'BetacslopeLSLH','BetacslopeRSLH','BetacslopeLSRH','BetacslopeRSRH',...
        'BetacampatresponseLSLH','BetacampatresponseRSLH','BetacampatresponseLSRH','BetacampatresponseRSRH'
 
        };
    
    IdlevelLong = {'Subject','group','bin','hand','RT','CPPonset','CPPslope','CPPampatresponse','CPP_T',...
        'Betacamp','BetacSlope'};
    
    cell2csv (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix150_bins.csv'],Idlevel) %create the csv file
    cell2csv (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix_R150_bins.csv'],IdlevelLong) %create the csv file
    
    subject_folder = {'AB99C','AG26C','AH05D','AL20D','AM23D','AZ60C','BB12D','BB47C','BD38D','BK48C',...
        'BL31D','CA30D','CC27D','CC76C','CD11D','CF80C','CG52C','CH78C','DN64C','DR69C',...
        'DR79C','DW45D','EB68D','ED62C','EF65C','EH43D','EV61C','FH003D','FO37D','GK67C',...
        'HT17D','HT72C','IA14D','IM19D','JC002D','JD35C','JD40C','JE04D','JE18D','JG10D',...
        'JG77C','JH63C','JK33D','JM25C','JM53D','JN21C','JP82C','JS24D','JS73C','JW71C',...
        'LB49C','LK32D','MB13D','MB16D','MD36C','MG01D','MK84C','MW70C','MY83C','NP44C',...
        'NV75C','OB50C','OH08D','OK34D','PD06D','RD15D','SH81C','SS85C','TB28C','TG51C',...
        'TK66C','VR86C','WH09D','WM22D','ZB41D'};
    cntLong=1;
    %%
    for c=1:3
        for s = 1:size(ssj{c},2)
            sT = ssj{c}(s);
            if CSD
                load(['Data/ERPs/' method '/CSD_' allnamefile{methodindx} '_' dyscAMRM{c} '_group_plots_2bins_' num2str(sT)],'CPP_side','CPPr_side','ERP_side','ERPr_side',...
                    'RTs','RT_bins','t','plot_chans','side_tags','tr','chanlocs',...
                    'hit','wrong','miss','allTs',...
                    'CPP_side_onsets','CPPr_slope','CPP_click','CPPr_peakTime',...
                    'RT_all','hit_rate','fAALL','fA_all','fA_count',...
                    'betac_side_onsets','betac_slope','Betac_click','betac_side','betarc_side',...
                    'betai_side_onsets','betai_slope','Betai_click','betai_side','betari_side',...
                    'bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side','STFT_timer','STFT_time')
            else
                load(['Data/ERPs/' method '/' dyscAMRM{c} 'group_plots_15_' num2str(sT)],'CPP_side','CPPr_side','ERP_side','ERPr_side',...
                    'N2c_side','N2i_side','RTs','t','plot_chans','side_tags','tr','ERP_hemi_side','ERPr_hemi_side','chanlocs',...
                    'hit','wrong','miss','allTs',...
                    'CPP_side_onsets','CPPr_slope','CPP_click','N2i_peak_amp_index_t','N2c_peak_amp_index_t','max_peak_N2i','max_peak_N2c','CPPr_peakTime',...
                    'RT_all','hit_rate','fAALL','fA_all','fA_count',...
                    'betac_side_onsets','betac_slope','Betac_click','betac_side','betarc_side',...
                    'betai_side_onsets','betai_slope','Betai_click','betai_side','betari_side',...
                    'bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side','STFT_timer','STFT_time',...
                    'alpha_side','alpha_base_side','alpha_asym_side','alpha_asym_avg_side','alpha_t','alphar_base_side','alphar_asym_side','alpha_tr',...
                    'preAlpha','preAlpha_L','preAlpha_R','preAlpha_asym','postAlpha','postAlpha_L','postAlpha_R','postAlpha_asym',...
                    'postAlphaR','postAlphaR_L','postAlphaR_R','postAlphaR_asym')
            end
            ID = subject_folder{sT};
            for hand=1:2
                for bins=1:no_of_bins
                    RT_stores{c,s,hand,bins} = RT_bins{1,hand,bins};
                    RT_all_bins(1,hand,bins) = nanmean([RT_bins{1,hand,bins}]);
                end
            end
            hit_index{c}(s,:) = 100*length([hit{:}])/(length([allTs{:}]));
            wrong_index{c}(s,:) = 100*length([wrong{:}])/(length([allTs{:}]));
            miss_index{c}(s,:) = 100*length([miss{:}])/(length([allTs{:}]));
            
            
            diffHitWrong{c}(s,:) = length([hit{:}])+length([wrong{:}])+length([miss{:}])-length([allTs{:}]);
            CPP_sides{c}(s,:,:,:)= squeeze(CPP_side);
            CPPr_sides{c}(s,:,:,:) = squeeze(CPPr_side);
            
            ERP_sides{c}(s,:,:,:,:) = squeeze(ERP_side);
            ERPr_sides{c}(s,:,:,:,:) = squeeze(ERPr_side);
            
            
            CPP_sides{c}(s,:,:,:)= squeeze(CPP_side);
            CPPr_sides{c}(s,:,:,:) = squeeze(CPPr_side);
            
            betac_sides{c}(s,:,:,:,:) = squeeze(betac_side);
            betarc_sides{c}(s,:,:,:,:) = squeeze(betarc_side);
            betai_sides{c}(s,:,:,:,:) = squeeze(betai_side);
            betari_sides{c}(s,:,:,:,:) = squeeze(betari_side);
            
            bERP_sides{c}(s,:,:,:,:) = squeeze(bERP_side_all);
            bERPr_sides{c}(s,:,:,:,:) = squeeze(bERPr_side_all);
            
            bERP_hemi_sides{c}(s,:,:,:,:) = squeeze(bERP_hemi_side);
            bERPr_hemi_sides{c}(s,:,:,:,:) = squeeze(bERP_r_hemi_side);
            
            % SSVEP            
            RT_sides{s,c}= [RTs{:}];
            RT_mean{c}(s) = mean([RT_sides{s,c}]);
            
            CPPclick_sides{c}(s,:,:) = squeeze(CPP_click);
            BetaSlope_sides{c}(s,:,:) = squeeze(betac_slope);
            
            
            temp = [fAALL{:}];
            fA_sides{s,c} = temp;
            fA_total{c}(s) = sum(fA_all(:));
            
            fA_aall{c}(s) = sum(fA_count(:));
            cplevelAll{c}(s) = mean(mean(CPP_click,2),3);
            cpslope{c}(s) = mean(mean(CPPr_slope,2),3);
            count = count+1;
            participant_level(count,:) = [sT,c,RT_all_bins(:)',CPP_side_onsets(:)',CPPr_slope(:)',...
                CPP_click(:)',CPPr_peakTime(:)',...
                betac_slope(:)',Betac_click(:)'
                ];
            
            for bin=1:no_of_bins
                for hand=1:2
                    participant_long(cntLong,:) = [sT,c,bin,hand,RT_all_bins(:,hand,bin)',CPP_side_onsets(:,hand,bin)',CPPr_slope(:,hand,bin)',...
                        CPP_click(:,hand,bin)',CPPr_peakTime(:,hand,bin)',...
                        Betac_click(:,hand,bin)',betac_slope(:,hand,bin)'
                        ];
                    cntLong=cntLong+1;
                end
            end
        end
    end
    %% Make participant level matrix for export into SPSS or R
    dlmwrite (savenameSPSS,participant_level,'-append')
    dlmwrite (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix_R150_bins.csv'],participant_long,'-append')
    
    %% plot them
    red=[1 0 0];blue=[0 0.45 0.74];black=[0 0 0];green=[0.47 0.67 0.19];
    colour = [red;blue;green];
    clear red blue yellow
    legend_tags = {'AM','RM','DD'};
    
    for c=1:3
        [~, indx{c}] = sort(RT_mean{c});
        fastindx{c} = indx{c}(1:ceil(length(indx{c})/2));
        slowindx{c} = indx{c}(ceil(length(indx{c})/2)+1:end);
    end
end
return
%% Behaviour
%% RT Behaviour
% make bar graphs
for bins=1:no_of_bins
    for c=1:3
        mean_RTs(c,bins) = mean([RT_stores{c,:,:,bins}]);
        count_RTs(c,bins) = length([RT_stores{c,:,:,bins}]);
        std_RTs(c,bins) = std([RT_stores{c,:,:,bins}])/sqrt(length(ssj{c}));
    end
end
hf = figure(9999)
hold on
cs=0;
for c=1:3
    hs1(c) = plot(1:no_of_bins,mean_RTs(c,:),'Linewidth',1.5);
    set(hs1(c),'Color',colour(c,:));
    hs2(c) = errorbar(1:no_of_bins,mean_RTs(c,:),std_RTs(c,:),'.');
    set(hs2(c),'Color',colour(c,:));
    hs2(c).LineWidth = 2;
end
ylabel({'Reaction Time (ms)'},'FontSize',16,'fontweight','bold')
axis([0 6 300 1200]);
xlabel({'Bins'},'FontSize',16,'fontweight','bold')
set(gca,'Fontsize',16,'fontweight','bold')
yticks(300:300:1200);
xticks(1:1:5);
legend(hs1,legend_tags,'FontSize',16,'Location','eastoutside','fontweight','normal')
saveas(hf,['Figures\NicoleGraphs\RT_bar_graphs_bins' method '.png']);
%%
hfig = figure(4)
for c=1:3
    %Beta scal topo
    t1 = -100; t2 = 100;
    plot_mean = squeeze(nanmean(nanmean(nanmean(nanmean(bERPr_hemi_sides{c}(:,:,find(STFT_timer>=t1 & STFT_timer<t2),:,:),1),3),4),5));

    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
          [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);    
end
suptitle('BetaPlots')
saveas(hfig,['Figures/NicoleGraphs/betatopoplot.jpg'])
%% Plot Beta x target side
for c=1:3
    Beta{c} = squeeze(nanmean(nanmean(betac_sides{c},1),3)); % time x side, no hand
    Betar{c} = squeeze(nanmean(nanmean(betarc_sides{c},1),3));
end

hfig=figure
clear h
for c=1:3
    
    for bin=1:no_of_bins
        
        h(c) = plot(STFT_timer,squeeze(Betar{c}(:,bin)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:200:100],'ylim',[-0.7,0.1],'ytick',[-0.7:0.2:0.1]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Beta'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
end
saveas(hfig,['Figures\NicoleGraphs\beta' method '.jpg'])
%%
% Grand average ERP
%chan x time x c
clear ERP_group_all
for c=1:3
ERP_group_all(:,:,c) = squeeze(mean(mean(mean(bERPr_hemi_sides{c}(:,:,:,:,:),1),4),5));
end

figure
plottopo(ERP_group_all(:,:,:),'chanlocs',chanlocs,'limits',[STFT_timer(1) STFT_timer(end) ...
    min(min(min(ERP_group_all(plot_chans,:,:))))  max(max(max(ERP_group_all(plot_chans,:,:))))], ...
    'title',['ERPr left vs right targets'],'legend',dyscAMRM,'showleg','on','ydir',1)
%%
for c=1:3
ERP_group_all(:,:,c) = squeeze(mean(mean(mean(ERPr_hemi_sides{c}(:,:,:,:,:),1),4),5));
end

figure
plottopo(ERP_group_all(:,:,:),'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
    min(min(min(ERP_group_all(plot_chans,:,:))))  max(max(max(ERP_group_all(plot_chans,:,:))))], ...
    'title',['ERPr left vs right targets'],'legend',dyscAMRM,'showleg','on','ydir',1)
%%
clear ERP_group_all
for c=1:3
ERP_group_all(:,:,c) = squeeze(mean(mean(mean(ERP_hemi_sides{c}(:,:,:,:,:),1),4),5));
end

figure
plottopo(ERP_group_all(:,:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    min(min(min(ERP_group_all(plot_chans,:,:))))  max(max(max(ERP_group_all(plot_chans,:,:))))], ...
    'title',['ERPr left vs right targets'],'legend',dyscAMRM,'showleg','on','ydir',1)
%% Plot CPP x target side
for c=1:3
    CPP{c} = squeeze(mean(CPP_sides{c},1)); % time x side, no hand
    CPPr{c} = squeeze(mean(CPPr_sides{c},1));
end

clear h
hfig= figure
clear h
for c=1:3
    for bin=1:no_of_bins
        
        h(c) = plot(t,squeeze(nanmean(CPP{c}(:,:,bin),2)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:200:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['CPP (Stimulus Locked)'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        
    end
end
saveas(hfig,['Figures\NicoleGraphs\CPP_' method 'all.jpg'])
%%
clear h
hfig= figure
for c=1:3
    for bin=1:no_of_bins
        h(c) = plot(tr,squeeze(nanmean(CPPr{c}(:,:,bin),2)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
    end
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['CPP (Response Locked)'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
end
legend(h,legend_tags, ...
    'FontSize',16,'Location','eastoutside');
saveas(hfig,['Figures\NicoleGraphs\CPPr_bins' method 'all.jpg'])
%%
% CPP amplitude across bins
for bins=1:no_of_bins
    for c=1:3
        mean_CPPamps(c,bins) = nanmean(nanmean([CPPclick_sides{c}(:,:,bins)],2));
        std_CPPamps(c,bins) = nanstd(nanmean([CPPclick_sides{c}(:,:,bins)],2))/sqrt(length(ssj{c}));
    end
end
hf = figure(9991)
hold on
cs=0;
for c=1:3
    hs1(c) = plot(1:no_of_bins,mean_CPPamps(c,:),'Linewidth',1.5);
    set(hs1(c),'Color',colour(c,:));
    hs2(c) = errorbar(1:no_of_bins,mean_CPPamps(c,:),std_CPPamps(c,:),'.');
    set(hs2(c),'Color',colour(c,:));
    hs2(c).LineWidth = 2;
end
set(gca,'Fontsize',16,'fontweight','bold')
ylabel('Resp CPP Amplitude (\muVolts)','FontSize',16,'fontweight','bold')
axis([0 6 0 70]);
xticks([1:1:5])
xlabel({'Bins'},'FontSize',16,'fontweight','bold')
yticks(0:20:100);
legend(hs1,legend_tags,'FontSize',16,'Location','eastoutside','fontweight','normal')
saveas(hf,['Figures\NicoleGraphs\CPPr_graphs_bins' method '.png']);

%%
% Beta slope across bins
for bins=1:no_of_bins
    for c=1:3
        mean_Betaamps(c,bins) = nanmean(nanmean([BetaSlope_sides{c}(:,:,bins)],2));
        std_Betaamps(c,bins) = nanstd(nanmean([BetaSlope_sides{c}(:,:,bins)],2))/sqrt(length(ssj{c}));
    end
end
hf = figure(9991)
hold on
cs=0;
for c=1:3
    hs1(c) = plot(1:no_of_bins,mean_Betaamps(c,:),'Linewidth',1.5);
    set(hs1(c),'Color',colour(c,:));
    hs2(c) = errorbar(1:no_of_bins,mean_Betaamps(c,:),std_Betaamps(c,:),'.');
    set(hs2(c),'Color',colour(c,:));
    hs2(c).LineWidth = 2;
end
ylabel('Resp Beta Slope (\muVolts/ms)','FontSize',16,'fontweight','bold')
axis([0 6 -0.003 0.001]);
xticks([1:1:5])
xlabel({'Bins'},'FontSize',16,'fontweight','bold')
legend(hs1,legend_tags,'Location','SouthWest')
set(gca,'Fontsize',16,'fontweight','bold')
yticks(-0.003:0.001:0.001);
legend(hs1,legend_tags,'FontSize',16,'Location','eastoutside','fontweight','normal')
saveas(hf,['Figures\NicoleGraphs\Beta_graphs_bins' method '.png']);
return

