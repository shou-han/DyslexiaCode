% establish the groups 
Methods={'HAPPE'};
allnamefile={'35'}; 
for methodindx=1
    method = Methods{methodindx};
    clearvars -except method Methods  methodindx allnamefile
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
    %ssj = 1;
    %ssj = [8 10 12:20 21 22:26 29];
    %addpath('function_programs/');
    %ssj = [1:5 8:26 28:29];
    CSD = 1;
    
    count = 0;
    Idlevel= {'Subject','group',...
        'RTLSLH','RTRSLH','RTLSRH','RTRSRH',...
        'AccuracyLSLH','AccuracyRSLH','AccuracyLSRH','AccuracyRSRH',...
        'falseAlarmLSLH','falseAlarmRSLH','falseAlarmLSRH','falseAlarmRSRH',...
        'CPPonsetLSLH','CPPonsetRSLH','CPPonsetLSRH','CPPonsetRSRH',...
        'CPPslopeLSLH','CPPslopeRSLH','CPPslopeLSRH','CPPslopeRSRH',...
        'CPPampatresponseLSLH','CPPampatresponseRSLH','CPPampatresponseLSRH','CPPampatresponseRSRH',...
        'CPPPeakTimeLSLH','CPPPeakTimeRSLH','CPPPeakTimeLSRH','CPPPeakTimeRSRH',...
        'N2ipeakTimeLSLH','N2ipeakTimeRSLH','N2ipeakTimeLSRH','N2ipeakTimeRSRH',...
        'N2cpeakTimeLSLH','N2cpeakTimeRSLH','N2cpeakTimeLSRH','N2cpeakTimeRSRH',...
        'N2ipeakLSLH','N2ipeakRSLH','N2ipeakLSRH','N2ipeakRSRH',...
        'N2cpeakLSLH','N2cpeakRSLH','N2cpeakLSRH','N2cpeakRSRH',...
        'BetaconsetLSLH','BetaconsetRSLH','BetaconsetLSRH','BetaconsetRSRH',...
        'BetacslopeLSLH','BetacslopeRSLH','BetacslopeLSRH','BetacslopeRSRH',...
        'BetacampatresponseLSLH','BetacampatresponseRSLH','BetacampatresponseLSRH','BetacampatresponseRSRH',...
        'BetaionsetLSLH','BetaionsetRSLH','BetaionsetLSRH','BetaionsetRSRH',...
        'BetaislopeLSLH','BetaislopeRSLH','BetaislopeLSRH','BetaislopeRSRH',...
        'BetaiampatresponseLSLH','BetaiampatresponseRSLH','BetaiampatresponseLSRH','BetaiampatresponseRSRH' , ...
        'preAlphaLSLH','preAlphaRSLH','preAlphaLSRH','preAlphaRSRH' , ...
        'preAlphaLLSLH','preAlphaLRSLH','preAlphaLLSRH','preAlphaLRSRH' , ...
        'preAlphaRLSLH','preAlphaRRSLH','preAlphaRLSRH','preAlphaRRSRH' , ...
        'preAlphaasymLSLH','preAlphaasymRSLH','preAlphaasymLSRH','preAlphaasymRSRH' , ...
        'postAlphaLSLH','postAlphaRSLH','postAlphaLSRH','postAlphaRSRH' , ...
        'postAlphaLLSLH','postAlphaLRSLH','postAlphaLLSRH','postAlphaLRSRH' , ...
        'postAlphaRLSLH','postAlphaRRSLH','postAlphaRLSRH','postAlphaRRSRH' , ...
        'postAlphaasymLSLH','postAlphaasymRSLH','postAlphaasymLSRH','postAlphaasymRSRH',...
        'postAlpharespLSLH','postAlpharespRSLH','postAlpharespLSRH','postAlpharespRSRH' , ...
        'postAlpharespLLSLH','postAlpharespLRSLH','postAlpharespLLSRH','postAlpharespLRSRH' , ...
        'postAlpharespRLSLH','postAlpharespRRSLH','postAlpharespRLSRH','postAlpharespRRSRH' , ...
        'postAlpharespasymLSLH','postAlpharespasymRSLH','postAlpharespasymLSRH','postAlpharespasymRSRH',...
        };
    IdlevelLong = {'Subject','group','hemi','hand','RT','Accuracy','falseAlarm','CPPonset','CPPslope','CPPampatresponse','CPP_T','N2i_T','N2c_T','N2iP','N2cP',...
        'Betacamp','BetacSlope','BetacOnset','Betaiamp','BetaiSlope','BetaiOnset',...
        'preAlpha','preAlphaL','preAlphaR','preAlphaasym','postAlpha','postAlphaL','postAlphaR','postAlphaasym',...
        'postAlpharesp','postAlphaLresp','postAlphaRresp','postAlphaasymresp'};
    
    cell2csv (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix150.csv'],Idlevel) %create the csv file
    cell2csv (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix_R150.csv'],IdlevelLong) %create the csv file
    
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
                load(['Data/ERPs/' method '/CSD_' allnamefile{methodindx} '_' dyscAMRM{c} '_group_plots_' num2str(sT)],'CPP_side','CPPr_side','ERP_side','ERPr_side',...
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
            
            hit_index{c}(s,:) = 100*length([hit{:}])/(length([allTs{:}]));
            wrong_index{c}(s,:) = 100*length([wrong{:}])/(length([allTs{:}]));
            miss_index{c}(s,:) = 100*length([miss{:}])/(length([allTs{:}]));
            
            
            diffHitWrong{c}(s,:) = length([hit{:}])+length([wrong{:}])+length([miss{:}])-length([allTs{:}]);
            CPP_sides{c}(s,:,:,:)= squeeze(CPP_side);
            CPPr_sides{c}(s,:,:,:) = squeeze(CPPr_side);
            N2c_sides{c}(s,:,:,:)=squeeze(N2c_side);
            N2i_sides{c}(s,:,:,:) =squeeze(N2i_side);
            
            ERP_sides{c}(s,:,:,:,:) = squeeze(ERP_side);
            ERPr_sides{c}(s,:,:,:,:) = squeeze(ERPr_side);
            
            ERP_hemi_sides{c}(s,:,:,:,:) = squeeze(ERP_hemi_side);
            ERPr_hemi_sides{c}(s,:,:,:,:) = squeeze(ERPr_hemi_side);
            
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
            
            alpha_sides{c}(s,:,:,:,:) = squeeze(alpha_side);
            alpha_base_sides{c}(s,:,:,:,:) = squeeze(alpha_base_side);
            alpha_asym_sides{c}(s,:,:,:,:) = squeeze(alpha_asym_side);
            alpha_asym_avg_sides{c}(s,:,:,:,:) = squeeze(alpha_asym_avg_side);
            alpha_r_sides{c}(s,:,:,:,:) = squeeze(alphar_base_side);
            alphar_asym_sides{c}(s,:,:,:,:) = squeeze(alphar_asym_side);
            
            
            RT_sides{s,c}= [RTs{:}];
            RT_mean{c}(s) = mean([RT_sides{s,c}]);
            
            CPPclick_sides{c}(s,:,:) = squeeze(CPP_click);
            
            
            temp = [fAALL{:}];
            fA_sides{s,c} = temp;
            fA_total{c}(s) = sum(fA_all(:));
            
            fA_aall{c}(s) = sum(fA_count(:));
            cplevelAll{c}(s) = mean(mean(CPP_click,2),3);
            alphalevelAll{c}(s,:,:) = preAlpha;
            alphalevelAllt{c}(s,:,:) = mean(mean(preAlpha,2),3);;
            cpslope{c}(s) = mean(mean(CPPr_slope,2),3);
            count = count+1;
            participant_level(count,:) = [sT,c,RT_all(:)',hit_rate(:)',fA_all(:)',CPP_side_onsets(:)',CPPr_slope(:)',...
                CPP_click(:)',CPPr_peakTime(:)',N2i_peak_amp_index_t(:)',N2c_peak_amp_index_t(:)',...
                max_peak_N2i(:)',max_peak_N2c(:)',...
                betac_side_onsets(:)',betac_slope(:)',Betac_click(:)',...
                betai_side_onsets(:)',betai_slope(:)',Betai_click(:)',...
                preAlpha(:)',preAlpha_L(:)',preAlpha_R(:)',preAlpha_asym(:)',...
                postAlpha(:)',postAlpha_L(:)',postAlpha_R(:)',postAlpha_asym(:)',...
                postAlphaR(:)',postAlphaR_L(:)',postAlphaR_R(:)',postAlphaR_asym(:)',...
                ];
            
            for side=1:2
                for hand=1:2
                    participant_long(cntLong,:) = [sT,c,side,hand,RT_all(:,side,hand)',hit_rate(:,side,hand)',fA_all(:,side,hand)',CPP_side_onsets(:,side,hand)',CPPr_slope(:,side,hand)',...
                        CPP_click(:,side,hand)',CPPr_peakTime(:,side,hand)',N2i_peak_amp_index_t(:,side,hand)',N2c_peak_amp_index_t(:,side,hand)',...
                        max_peak_N2i(:,side,hand)',max_peak_N2c(:,side,hand)',...
                        Betac_click(:,side,hand)',betac_slope(:,side,hand)',betac_side_onsets(:,side,hand)',...
                        Betai_click(:,side,hand)',betai_slope(:,side,hand)',betai_side_onsets(:,side,hand)',...
                        preAlpha(:,side,hand)',preAlpha_L(:,side,hand)',preAlpha_L(:,side,hand)',preAlpha_asym(:,side,hand)',...
                        postAlpha(:,side,hand)',postAlpha_L(:,side,hand)',postAlpha_R(:,side,hand)',postAlpha_asym(:,side,hand)',...
                        postAlphaR(:,side,hand)',postAlphaR_L(:,side,hand)',postAlphaR_R(:,side,hand)',postAlphaR_asym(:,side,hand)',...
                        ];
                    cntLong=cntLong+1;
                end
            end
        end
    end
    %% Make participant level matrix for export into SPSS or R
    dlmwrite (savenameSPSS,participant_level,'-append')
    dlmwrite (['../Analyses_Scripts_R/Stats/' method '_participant_level_matrix_R150.csv'],participant_long,'-append')
    
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
%%
for c=1:3
    hold all
    plot3(cpslope{c},cplevelAll{c},fA_aall{c},'+');
end
%%
for c=1:3
    % Quick check of RT index
    [~,p,~,stats] = ttest(RT_indexs{c});
    disp(['RT index to zero: t = ' num2str(stats.tstat) ', p = ' num2str(p)])
    
    figure
    plot(zscore(RT_indexs{c}))
end
%% Behaviour
%% RT Behaviour
clear hs hb
    for c = 1:3
        temp  =[];
        for s = 1:length(ssj{c})
            temp = [temp [RT_sides{s,c}]];
        end
        meanRTs{c} = temp;
        clear temp 
    end
hf = figure(9998)
hold on
cs=0;
for c = 1:3
    cs=cs+1;
    clear de_temp
    de_temp = [RT_sides{:,c}];
      set(gca,'FontSize',16,'xlim',[300,1600],'xtick',[300:200:1600],'ylim',[0,600],'ytick',[0:200:600]);%,'ylim',[-1.5,0.5]);
    h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
    hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',colour(c,:),'Linewidth',1.5);
    hb(cs) = line([mean(meanRTs{c}) mean(meanRTs{c})],[0 450],'Color',colour(c,:),'Linewidth',1.5,'LineStyle','--');    
end
ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
legend(hs,legend_tags,'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\NicoleGraphs\RT_all_bins' method '.png']);
% make bar graphs
for c=1:3
    mean_RTs(c) = mean([RT_sides{:,c}]);
    std_RTs(c) = std([RT_sides{:,c}])/sqrt(length(ssj{c}));
end
hf = figure(9999)
hold on
cs=0;
for c=1:3
    hs1(c) = bar(c,mean_RTs(c),'Linewidth',1.5);
    set(hs1(c),'FaceColor',colour(c,:));
    hs2(c) = errorbar(c,mean_RTs(c),std_RTs(c)'.');
    hs2(c).Color = [0 0 0];
    hs2(c).LineWidth = 2;33
end
ylabel({'Response Time (ms)'},'FontSize',16,'fontweight','bold')
axis([0 4 300 900]);33
xticks([1 2 3]);xticklabels(legend_tags);
set(gca,'Fontsize',20)
yticks(300:300:1200);
legend(hs,legend_tags,'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\NicoleGraphs\RT_bar_graphs' method '.png']);

%% False Alarm
% make bar graphs
for c=1:3
    mean_FA(c) = mean(fA_total{c});
    std_FA(c) = std(fA_total{c})/sqrt(length(ssj{c}));
end
hf = figure(9999)
hold on
cs=0;
for c=1:3
    hs1(c) = bar(c,mean_FA(c),'Linewidth',1.5);
    set(hs1(c),'FaceColor',colour(c,:));
    hs2(c) = errorbar(c,mean_FA(c),std_FA(c)'.');
    hs2(c).Color = [0 0 0];
    hs2(c).LineWidth = 2;
end
ylabel({'False Alarm'},'FontSize',16,'fontweight','bold')
axis([0 4 0 20]);
xticks([1 2 3]);xticklabels(legend_tags);
set(gca,'Fontsize',20)
yticks(0:5:20);
legend(hs1,legend_tags,'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\NicoleGraphs\FA_bar_graphs' method '.png']);
%%
clear alpha_asym_group_side
%RM-AM
for c=1:3
alpha_asym_group_side(:,:,c) = squeeze(mean(mean(alpha_sides{c}(:,:,:,:),1),4)); % chan x time x side
end
for side=1
    figure
plottopo(alpha_asym_group_side(:,:,:),'chanlocs',chanlocs,'limits',[-1000 1543 ...
    min(min(min(alpha_asym_group_side(1:64,:,:))))  max(max(max(alpha_asym_group_side(1:64,:,:))))], ...
    'title',['Alpha asymmetry left vs right targets'],'legend',side_tags,'showleg','on','ydir',1)
end

%% topoplot alphas
for c=1:3
clear plot_mean
t1 = -300; t2 = 100;
hfAlpha = figure(92)
subplot(1,3,c)
plot_mean = squeeze(mean(mean(mean(alpha_sides{c}(:,:,find(alpha_t>=t1 & alpha_t<=t2),:),1),3),4));
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean) max(plot_mean)], ...
    'electrodes','off','plotchans',plot_chans);
title(dyscAMRM{c})
end
saveas(hfAlpha,['Figures\alpha' method '.png']);
%%
clear alpha_base_group_side
for c=1:3
t1 = -100; t2 = 100;
alpha_base_group_side = squeeze(mean(alphar_asym_sides{c},1));
for side = 1
    plot_mean = squeeze(mean(mean(alpha_base_group_side(:,find(alpha_tr>=t1 & alpha_tr<=t2),:),2),3));
    halpha = figure
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [-0.07 0.01], ...
        'electrodes','numbers','plotchans',plot_chans);    
    title(['PostTarget Alpha Resp -100ms before to 100ms after Targets'],'FontSize',12);
end
suptitle([dyscAMRM{c}]);
saveas(halpha,['Figures\alpha' num2str(c) '.png']);
end

%%
figure
for c=1:3
t1 = 500; t2 = 700;
alpha_asym_group_side = squeeze(mean(alpha_asym_sides{c},1));
for side = 1
    plot_mean = squeeze(mean(mean(alpha_asym_group_side(:,find(alpha_t>=t1 & alpha_t<=t2),:),2),3));
    figure
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [-0.07, ...
        0], ...
        'electrodes','numbers','plotchans',plot_chans);    
    colorbar('FontSize',12)
    title(['PostTarget Alpha Asym 150ms to 700ms Targets'],'FontSize',12);
end
suptitle([dyscAMRM{c}]);
end

%%
clear alpha_base_group_side
figure
for c=1:3
t1 = 300; t2 = 700;
alpha_base_group_side = squeeze(mean(alpha_base_sides{c},1));
for side = 1
    plot_mean = squeeze(mean(mean(alpha_base_group_side(:,find(alpha_t>=t1 & alpha_t<=t2),:),2),3));
    halpha=figure
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(squeeze(mean(mean(alpha_base_group_side(:,find(alpha_t>=t1 & alpha_t<=t2),:),3),2))) ...
        max(squeeze(mean(mean(alpha_base_group_side(:,find(alpha_t>=t1 & alpha_t<=t2),:),3),2)))], ...
        'electrodes','off','plotchans',plot_chans);    
    colorbar('FontSize',12)
    title(['PostTarget Alpha 150ms to 700ms Targets'],'FontSize',12);
end
suptitle([dyscAMRM{c}]);
saveas(halpha,['Figures\alphapost' num2str(c) '.png']);
end

%% plot the three alpha waves
clear pre_alpha_hemi hss hfig
 ch_alpha{1} = [61 62 63]; ch_alpha{2} =[61 62 63] ;
for c=1:3
    for hemi=1:2
        pre_alpha_hemi(c,hemi,:)=squeeze(mean(mean(mean(mean(alpha_sides{c}(:,ch_alpha{1},:,:,:),1),4),2),5));
    end
end
for hemi=1
    hfig(hemi) = figure(hemi)
    hold on
    for c=1:3
        allData = squeeze(mean(mean(mean(alpha_sides{c}(:,ch_alpha{1},:,:,:),5),4),2));
        meanallData = mean(allData);
        stdallData = std(allData)/sqrt(size(allData,1));
        hE(c)=shadedErrorBar(alpha_t,meanallData,stdallData,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});hold on;    
        hss(c)=plot(alpha_t,squeeze(mean(pre_alpha_hemi(c,hemi,:),2)),'Color',colour(c,:),'LineWidth',3,'LineStyle','-')
        set(gca,'FontSize',16,'xlim',[-500,1200],'xtick',[-500, 0:400:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Pre-Target Alpha'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
        legend(hss,legend_tags, ...
            'FontSize',16,'Location','eastoutside');    
saveas(hfig(hemi),['Figures\NicoleGraphs\AlphaReview' method '.jpg'])
end
%% plot the three alpha waves Reviewer
clear pre_alpha_hemi hss hfig
 ch_alpha{1} = [61 62 63]; ch_alpha{2} =[61 62 63] ;
for hemi=1
    hfig(hemi) = figure(hemi)
    hold on
    for c=1:3
        Alphadata = squeeze(mean(mean(mean(alpha_sides{c}(:,ch_alpha{2},:,:,:),2),4),5));
        meanAlpha = mean(Alphadata);
        stdAlpha = std(Alphadata,1)/sqrt(size(Alphadata,1));
        hold on
        hE(c) = shadedErrorBar(alpha_t,meanAlpha,stdAlpha,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});
        
        hss(c)=plot(alpha_t,meanAlpha,'Color',colour(c,:),'LineWidth',3,'LineStyle','-')
        set(gca,'FontSize',16,'xlim',[-500,1200],'xtick',[-500, 0:400:1200]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Pre-Target Alpha'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
        legend(hss,legend_tags, ...
            'FontSize',16,'Location','eastoutside');    
saveas(hfig(hemi),['Figures\NicoleGraphs\Alpha' method '_reviewer.jpg'])
end
%%
clear pre_alpha_hemi
 ch_alpha{1} = [59 64]; ch_alpha{2} =[59 64] ;
for c=1:3
    for hemi=1
        pre_alpha_hemi(c,hemi,:)=squeeze(mean(mean(mean(mean(alphar_asym_sides{c}(:,ch_alpha{hemi},:,:,:),1),5),4),2));
    end
end

for hemi=1
    hfig(hemi) = figure(hemi)
    hold on
    for c=1:3
        hss(c)=plot(alpha_tr,squeeze(pre_alpha_hemi(c,hemi,:,:)),'Color',colour(c,:),'LineWidth',3,'LineStyle','-')
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:200:100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Alpha Asymmetry'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
        legend(hss,legend_tags, ...
            'FontSize',16,'Location','NorthWest');    
saveas(hfig(hemi),['Figures\Alpha_asym.jpg'])
end
%%
clear pre_alpha_hemi
clear hss
 ch_alpha{1} = [61]; ch_alpha{2} =[61] ;
for c=1:3
    for hemi=1
        pre_alpha_hemi(c,hemi,:)=squeeze(mean(mean(mean(alpha_r_sides{c}(:,ch_alpha{hemi},:,hemi,:),1),5),2));
    end
end

for hemi=1
    hfig(hemi) = figure(hemi)
    hold on
    for c=1:3
        hss(c)=plot(alpha_tr,squeeze(mean(pre_alpha_hemi(c,hemi,:),2)),'Color',colour(c,:),'LineWidth',3,'LineStyle','-')
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:100:100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Post-Target Alpha'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
        legend(hss,legend_tags, ...
            'FontSize',16,'Location','eastoutside');    
saveas(hfig(hemi),['Figures\NicoleGraphs\Alpha_resp.jpg'])
end
%% Make some scalp plots for R:
%CPP scalp topo
for c=1:3
    t1 = -50; t2 = 50;
    plot_mean = squeeze(mean(mean(mean(mean(ERPr_hemi_sides{c}(:,:,find(tr>=t1 & tr<t2),:,:),1),3),4),5));
    scap1 = figure(1)
    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);
    title(legend_tags{c})
    
    %N2 scalp topo
    t1 = 280; t2 = 320;
    plot_mean = squeeze(mean(mean(mean(mean(ERP_hemi_sides{c}(:,:,find(t>=t1 & t<t2),:,:),1),3),4),5));
    scap2 = figure(2)
    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);
     title(legend_tags{c})
    %N2 scalp topo
    t1 = 280; t2 = 320;
    plot_mean = squeeze(mean(mean(mean(mean(ERP_hemi_sides{c}(:,:,find(t>=t1 & t<t2),:,:),1),3),4),5));
    scap3 = figure(3)
    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);
     title(legend_tags{c})
end
figure(1);suptitle('CPP');saveas(scap1, ['Figures/NicoleGraphs/cpptopoplot' method '.jpg'])
figure(2);suptitle('N2c');saveas(scap2, ['Figures/NicoleGraphs/N2ctopoplot' method '.jpg'])
figure(3);suptitle('N2i');saveas(scap3, ['Figures/NicoleGraphs/N2itopoplot' method '.jpg'])
%%
hfig = figure(4)
for c=1:3
    %Beta scal topo
    t1 = -100; t2 = 100;
    plot_mean = squeeze(mean(mean(mean(mean(bERPr_hemi_sides{c}(:,:,find(STFT_timer>=t1 & STFT_timer<t2),:,:),1),3),4),5));

    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
          [min(plot_mean) max(plot_mean)], ...
        'electrodes','off','plotchans',plot_chans);    
end
suptitle('BetaPlots')
saveas(hfig,['Figures/NicoleGraphs/betatopoplot.jpg'])
%% Plot Beta x target side
for c=1:3
    Beta{c} = squeeze(mean(mean(betac_sides{c},1),3)); % time x side, no hand
    Betar{c} = squeeze(mean(mean(betarc_sides{c},1),3));
end

clear h

for side=1
    hfig(side)=figure
    clear h
    for c=1:3
        h(c) = plot(STFT_timer,squeeze(mean(mean(Betar{c},2),4)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:200:100],'ylim',[-0.7,0.1],'ytick',[-0.7:0.2:0.1]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Beta'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        legend(h,legend_tags, ...
            'FontSize',16,'Location','eastoutside');
    end
       saveas(hfig(side),['Figures\NicoleGraphs\beta' method '.jpg'])
end
%% Plot Beta x target side
for c=1:3
    Beta{c} = squeeze(mean(mean(betac_sides{c},1),3)); % time x side, no hand
    Betar{c} = squeeze(mean(mean(betarc_sides{c},1),3));
end

clear h

for side=1
    hfig(side)=figure
    clear h
    for c=1:3
        betaData = squeeze(mean(mean(betarc_sides{c}(:,:,:,:),3),4));
        meanBeta = mean(betaData);
        stdBeta = std(betaData,1)/sqrt(size(betaData,1));
        hold on
        hE(c) = shadedErrorBar(STFT_timer,meanBeta,stdBeta,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});        
        h(c) = plot(STFT_timer,squeeze(mean(mean(Betar{c},2),4)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:200:100],'ylim',[-0.7,0.1],'ytick',[-0.7:0.2:0.1]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['Beta'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        legend(h,legend_tags, ...
            'FontSize',16,'Location','eastoutside');
    end
       saveas(hfig(side),['Figures\NicoleGraphs\betaReviewer' method '.jpg'])
end
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

for hemi=1
    hfig(hemi)= figure
    clear h
    for c=1:3
        %h(c) = plot(t,mean(mean(squeeze(CPP{c}(:,:,:)),2),3),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        CPPdata = squeeze(mean(mean(CPP_sides{c}(:,:,:),3),4));
        meanCPP = mean(CPPdata);
        stdCPP = std(CPPdata,1)/sqrt(size(CPPdata,1));
                hold on
        hE(c) = shadedErrorBar(t,meanCPP,stdCPP,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});
        h(c) = plot(t,mean(mean(squeeze(CPP{c}(:,:,:)),2),3),'LineWidth',3,'LineStyle','-','Color',colour(c,:))
        set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[0:400:1200],'ylim',[-4,10],'ytick',[-4:2:10]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['CPP (Stimulus Locked)'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    end
    legend(h,legend_tags, ...
    'FontSize',16,'Location','eastoutside');
    saveas(hfig(hemi),['Figures\NicoleGraphs\CPP_' method 'all.jpg'])
end
%%
for hemi=1
    hfig(hemi) = figure
    clear h
    for c = 1:3
        CPPrdata = squeeze(mean(mean(CPPr_sides{c}(:,:,:),3),4));
        meanCPPr = mean(CPPrdata);
        stdCPPr = std(CPPrdata,1)/sqrt(size(CPPrdata,1));
        hold on
        hE2(c) = shadedErrorBar(tr,meanCPPr,stdCPPr,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});
        h(c) = plot(tr,squeeze(mean(mean(CPPr{c}(:,:,:),2),3)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));
    end
    set(gca,'FontSize',16,'xlim',[-500,100],'ylim',[-4,10],'ytick',[-4:2:10]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['CPP (Response Locked)'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,legend_tags, ...
        'FontSize',16,'Location','eastoutside');
    saveas(hfig(hemi),['Figures\NicoleGraphs\CPPr_' method 'all.jpg'])
end




%% Plot N2c x target side
clear h
for c=1:3
    N2c{c} = squeeze(mean(ERP_hemi_sides{c}(:,23,:,:,:),1)); % time x side
end

for hand=1
    hfig(hand)= figure
    clear h
    for c=1:3
        N2cdata = squeeze(mean(mean(ERP_hemi_sides{c}(:,23,:,:,:),4),5));
        meanN2c = mean(N2cdata);
        stdN2c = std(N2cdata,1)/sqrt(size(N2cdata,1));
                hold on
        hE(c) = shadedErrorBar(t,meanN2c,stdN2c,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});

        h(c) = plot(t,squeeze(mean(mean(N2c{c}(:,:,:),2),3)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));
        hold on
    end
    set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:300:1200],'ylim',[-6,8],'ytick',[-6:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['N2c'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,legend_tags, ...
        'FontSize',16,'Location','eastoutside');
               saveas(hfig(hand),['Figures\NicoleGraphs\N2c' method '.jpg'])
end


%% Plot N2i x target side
for c=1:3
    N2i{c} = squeeze(mean(ERP_hemi_sides{c}(:,27,:,:,:),1)); % time x side
end
for hand=1
    hfig(hand)= figure
    clear h
    for c=1:3
       N2idata = squeeze(mean(mean(ERP_hemi_sides{c}(:,27,:,:,:),4),5));
        meanN2i = mean(N2idata);
        stdN2i = std(N2idata,1)/sqrt(size(N2idata,1));
                hold on
        hE(c) = shadedErrorBar(t,meanN2i,stdN2i,'lineprops',{'Color',colour(c,:),'LineWidth',3,'LineStyle','-'});

        h(c) = plot(t,squeeze(mean(mean(N2i{c}(:,:,:),2),3)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));
        hold on
    end
    set(gca,'FontSize',16,'xlim',[-100,1200],'xtick',[-100,0:300:1200],'ylim',[-4,3],'ytick',[-4:1:3])%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['N2i'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,legend_tags, ...
        'FontSize',16,'Location','eastoutside');
               saveas(hfig(hand),['Figures\NicoleGraphs\N2i' method '.jpg'])
end

return

