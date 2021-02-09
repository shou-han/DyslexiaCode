% establish the groups 
clear all
clc
close all
set(0,'DefaultFigureVisible','on')
STFT_timer=[];
addpath(genpath('../CSDtoolbox/'));
addpath(genpath('../eeglab13_6_5b/'));
addpath('functions/');
eeglab
hzusedTag = 1;
HzAll = [21,42];
HzAllName = {'21','42'};
Hzused = HzAll(hzusedTag);
%4,63,65, 19, 
ssj{1} = [1 2 10 16 19 21 24 27 30 32 36 37 41 42 44 46  51 60 62 67 68];%AM 47 58
ssj{2} = [6 8  17 18 20 25 49 50 55 57 59 61 69 71 72];%RM 14
ssj{3} = [3 5 7 9 11 12 13 15 22 23 26 28 29 31 33 34  38  40 43 45 48 52  54  53 56 35 39 66 64  73 74 75];%Dysl 

savenameSPSS = ['../Analyses_Scripts_R/Stats/participant_level_matrixSSVEP_' num2str(Hzused) '.csv']; 
dyscAMRM = {'AM_control','RM_control','Dyslexic'};
CSD = 1;

count = 0;
Idlevel= {'Subject','group',...
        'RTLSLH','RTRSLH','RTLSRH','RTRSRH',...
        'AccuracyLSLH','AccuracyRSLH','AccuracyLSRH','AccuracyRSRH',...          
        'SSVEPslopeLSLH','SSVEPslopeRSLH','SSVEPslopeLSRH','SSVEPslopeRSRH',...
        'SSVEPampatresponseLSLH','SSVEPampatresponseRSLH','SSVEPampatresponseLSRH','SSVEPampatresponseRSRH'
        };        
IdlevelLong = {'Subject','group','hemi','hand','RT','Accuracy','falseAlarm','CPPonset','CPPslope','CPPampatresponse','CPP_T','N2i_T','N2c_T','N2iP','N2cP',...
    'Betacamp','BetacSlope','BetacOnset','Betaiamp','BetaiSlope','BetaiOnset',...
    'preAlpha','preAlphaL','preAlphaR','preAlphaasym','postAlpha','postAlphaL','postAlphaR','postAlphaasym'};

cell2csv (savenameSPSS,Idlevel) %create the csv file
cell2csv (['../Analyses_Scripts_R/Stats/participant_level_matrix_R' num2str(Hzused) '.csv'],IdlevelLong) %create the csv file

subject_folder = {'AB99C','AG26C','AH05D','AL20D','AM23D','AZ60C','BB12D','BB47C','BD38D','BK48C',...
                  'BL31D','CA30D','CC27D','CC76C','CD11D','CF80C','CG52C','CH78C','DN64C','DR69C',...
                  'DR79C','DW45D','EB68D','ED62C','EF65C','EH43D','EV61C','FH003D','FO37D','GK67C',...
                  'HT17D','HT72C','IA14D','IM19D','JC002D','JD35C','JD40C','JE04D','JE18D','JG10D',...
                  'JG77C','JH63C','JK33D','JM25C','JM53D','JN21C','JP82C','JS24D','JS73C','JW71C',...
                  'LB49C','LK32D','MB13D','MB16D','MD36C','MG01D','MK84C','MW70C','MY83C','NP44C',...
                  'NV75C','OB50C','OH08D','OK34D','PD06D','RD15D','SH81C','SS85C','TB28C','TG51C',...
                  'TK66C','VR86C','WH09D','WM22D','ZB41D'};
cntLong=1;
for c=1:3
    for s = 1:size(ssj{c},2)
        sT = ssj{c}(s);
        if CSD
            load(['Data/ERPs/SSVEP_' num2str(Hzused) '_' dyscAMRM{c} 'group_plots_CSD' num2str(sT)],'CPP_side','CPPr_side','ERP_side','ERPr_side',...
                'N2c_side','N2i_side','RTs','t','plot_chans','side_tags','tr','ERP_hemi_side','ERPr_hemi_side','chanlocs',...
                'hit','wrong','miss','allTs',...
                'CPP_side_onsets','CPPr_slope','CPP_click','N2i_peak_amp_index_t','N2c_peak_amp_index_t','max_peak_N2i','max_peak_N2c','CPPr_peakTime',...
                'RT_all','hit_rate','fAALL','fA_all','fA_count',...
                'betac_side_onsets','betac_slope','Betac_click','betac_side','betarc_side',...
                'betaasym_side','betarasym_side','betai_side_onsets','betai_slope','Betai_click','betai_side','betari_side',...
                'bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side','STFT_timer','STFT_time',...
                'alpha_side','alpha_base_side','alpha_asym_side','alpha_asym_avg_side','alpha_t',...
                'preAlpha','preAlpha_L','preAlpha_R','preAlpha_asym','postAlpha','postAlpha_L','postAlpha_R','postAlpha_asym')
        
        else
            load(['Data/ERPs/' dyscAMRM{c} 'group_plots_SSVEP_60_20_6degree' num2str(sT)],'CPP_side','CPPr_side','ERP_side','ERPr_side',...
                'N2c_side','N2i_side','RTs','t','plot_chans','side_tags','tr','ERP_hemi_side','ERPr_hemi_side','chanlocs',...
                'hit','wrong','miss','allTs',...
                'CPP_side_onsets','CPPr_slope','CPP_click','N2i_peak_amp_index_t','N2c_peak_amp_index_t','max_peak_N2i','max_peak_N2c','CPPr_peakTime',...
                'RT_all','hit_rate','fAALL','fA_all','fA_count',...
                'betac_side_onsets','betac_slope','Betac_click','betac_side','betarc_side',...
                'betaasym_side','betarasym_side','betai_side_onsets','betai_slope','Betai_click','betai_side','betari_side',...
                'beta_asym','bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side','STFT_timer','STFT_time',...
                'alpha_side','alpha_base_side','alpha_asym_side','alpha_asym_avg_side','alpha_t',...
                'preAlpha','preAlpha_L','preAlpha_R','preAlpha_asym','postAlpha','postAlpha_L','postAlpha_R','postAlpha_asym')         
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
        bERPr_hemi_asym{c}(s,:,:,:,:) = squeeze(betarasym_side);
        alpha_sides{c}(s,:,:,:,:) = squeeze(alpha_side);
        alpha_base_sides{c}(s,:,:,:,:) = squeeze(alpha_base_side);
        alpha_asym_sides{c}(s,:,:,:,:) = squeeze(alpha_asym_side);
        alpha_asym_avg_sides{c}(s,:,:,:,:) = squeeze(alpha_asym_avg_side);

        
        RT_sides{s,c}= [RTs{:}];
        RT_mean{c}(s) = mean([RT_sides{s,c}]);
        
        CPPclick_sides{c}(s,:,:) = squeeze(CPP_click);
        
        
        temp = [fAALL{:}];
        fA_sides{s,c} = temp;
        fA_total{c}(s) = sum(fA_all(:));
        
        fA_aall{c}(s) = sum(fA_count(:));
        
        betaClick{c}(s,:,:) = Betac_click(:);
        betaiClick{c}(s,:,:) = Betai_click(:);
        cplevelAll{c}(s) = mean(mean(CPP_click,2),3);
        cpslope{c}(s) = mean(mean(CPPr_slope,2),3);
        count = count+1;
        participant_level(count,:) = [sT,c,RT_all(:)',hit_rate(:)',betac_slope(:)',Betac_click(:)'
            ];
        
        for side=1:2
            for hand=1:2
                participant_long(cntLong,:) = [sT,c,side,hand,RT_all(:,side,hand)',hit_rate(:,side,hand)',fA_all(:,side,hand)',CPP_side_onsets(:,side,hand)',CPPr_slope(:,side,hand)',...
                    CPP_click(:,side,hand)',CPPr_peakTime(:,side,hand)',N2i_peak_amp_index_t(:,side,hand)',N2c_peak_amp_index_t(:,side,hand)',...
                    max_peak_N2i(:,side,hand)',max_peak_N2c(:,side,hand)',...
                    Betac_click(:,side,hand)',betac_slope(:,side,hand)',betac_side_onsets(:,side,hand)',...
                    Betai_click(:,side,hand)',betai_slope(:,side,hand)',betai_side_onsets(:,side,hand)',...
                    preAlpha(:,side,hand)',preAlpha_L(:,side,hand)',preAlpha_L(:,side,hand)',preAlpha_asym(:,side,hand)',...
                    postAlpha(:,side,hand)',postAlpha_L(:,side,hand)',postAlpha_R(:,side,hand)',postAlpha_asym(:,side,hand)'];
                cntLong=cntLong+1;
            end
        end
    end
end    
%% Make participant level matrix for export into SPSS or R
dlmwrite (savenameSPSS,participant_level,'-append')
dlmwrite (['../Analyses_Scripts_R/Stats/participant_level_matrix_R' num2str(Hzused) '.csv'],participant_long,'-append')

%% plot them
red=[1 0 0];blue=[0 0.45 0.74];black=[0 0 0];green=[0.47 0.67 0.19];
colour = [red;blue;green];
clear red blue yellow
legend_tags = {'AM','RM','DD'};
leftright = {'left','right'};
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
        for s = 1:length(ssj)
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
      set(gca,'FontSize',16,'xlim',[300,1600],'xtick',[300:200:1600],'ylim',[0,800],'ytick',[0:200:800]);%,'ylim',[-1.5,0.5]);
    h = histc(de_temp, [min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)]);
    hs(cs) =plot([min(de_temp):(max(de_temp)-min(de_temp))/20:max(de_temp)],h,'Color',colour(c,:),'Linewidth',1.5);
    hb(cs) = line([mean(meanRTs{c}) mean(meanRTs{c})],[0 450],'Color',colour(c,:),'Linewidth',1.5,'LineStyle','--');    
end
ylabel({'Frequency across';' all Participants'},'FontSize',16,'fontweight','bold')
xlabel('Response Time (ms)','FontSize',16,'fontweight','bold')
legend(hs,legend_tags,'FontSize',12,'Location','NorthEast','fontweight','bold')
saveas(hf,['Figures\RT_all_bins.png']);
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
    set(hs1(c),'FaceColor',colour(c,:));33
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
saveas(hf,['Figures\RT_bar_graphs.png']);

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
saveas(hf,['Figures\FA_bar_graphs.png']);
%%
hfig = figure(4)
for c=1:3
    %Beta scal topo
    t1 = -100; t2 = 100;
    plot_mean = squeeze(mean(mean(mean(mean(bERPr_hemi_sides{c}(:,:,find(STFT_timer>=t1 & STFT_timer<t2),:,:),1),3),4),5));

    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(min(plot_mean)),max(max(plot_mean))], ...
        'electrodes','off','plotchans',plot_chans);    
end
colorbar
suptitle(['SSVEP '  HzAllName{hzusedTag} ' Scalp Plot'])
saveas(hfig,['Figures\NicoleGraphs\SSVEP_' HzAllName{hzusedTag} '_topo.jpg'])

%%
hfig = figure(4)
for c=1:3
    %Beta scal topo
    t1 = -100; t2 = 100;
    plot_mean = squeeze(mean(mean(mean(mean(bERPr_hemi_asym{c}(:,:,find(STFT_timer>=t1 & STFT_timer<t2),:,:),1),3),4),5));

    subplot(1,3,c)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(min(plot_mean)),max(max(plot_mean))], ...
        'electrodes','numbers','plotchans',plot_chans);    
end
colorbar
suptitle(['SSVEP '  HzAllName{hzusedTag} ' Scalp Plot'])
saveas(hfig,['Figures\NicoleGraphs\SSVEP_' HzAllName{hzusedTag} '_topoasym.jpg'])
%%
% Grand average ERP
%chan x time x c
for c=1:3
ERP_group_all(:,:,c) = squeeze(mean(mean(mean(bERPr_hemi_asym{c}(:,:,:,:,:),1),4),5));
end

figure
plottopo(ERP_group_all(:,:,:),'chanlocs',chanlocs,'limits',[STFT_time(1) STFT_time(end) ...
    min(min(min(ERP_group_all(plot_chans,:,:))))  max(max(max(ERP_group_all(plot_chans,:,:))))], ...
    'title',['ERPr left vs right targets'],'legend',dyscAMRM,'showleg','on','ydir',1)


%% Plot SSVEP x target side
for c=1:3
    ssvep{c} = squeeze(mean(betai_sides{c},1));
    ssvepr{c} = squeeze(mean(betari_sides{c},1));
end

clear h

for side=1
    hfig(side)=figure
    clear h
    for c=1:3
        h(c) = plot(STFT_time,squeeze(mean(mean(ssvep{c}(:,:,:),3),2)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        
        set(gca,'FontSize',12,'xlim',[-200,1500],'xtick',[-200,0:500:1500]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['SSVEP ' HzAllName{hzusedTag} 'Hz'])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        legend(h,legend_tags, ...
            'FontSize',16,'Location','eastoutside');
    end
       saveas(hfig(side),['Figures\NicoleGraphs\SSVEP' HzAllName{hzusedTag} '.jpg'])
end


for side=1
    hfig(side)=figure
    clear h
    for c=1:3
        h(c) = plot(STFT_timer,squeeze(mean(mean(ssvepr{c}(:,:,:),3),2)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
        
        set(gca,'FontSize',16,'xlim',[-500,100],'xtick',[-500:200:100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
        ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
        xlabel('Time (ms)','FontName','Arial','FontSize',16)
        title(['SSVEP ' HzAllName{hzusedTag} 'Hz (resp) '])
        line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
        line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        legend(h,legend_tags, ...
            'FontSize',16,'Location','eastoutside');
    end
       saveas(hfig(side),['Figures\NicoleGraphs\SSVEPReview' HzAllName{hzusedTag} '.jpg'])
end
return

