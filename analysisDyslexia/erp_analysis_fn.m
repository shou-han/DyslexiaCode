function erp_analysis_fn(number)
%number=6;
subj=1:75;
LPF=[35];
HPF=[5000];%5000 is HAPPE 
[a b c] = ndgrid(subj, LPF, HPF);
allcombinations = [a(:) b(:) c(:)];
conditions = allcombinations(number,:);
subj = conditions(1); LPFused=conditions(2);HPFused = conditions(3);jobID=number;
%%
clearvars -except subj LPFused HPFused pupil pupilr jobID
%% Use Current Source Density transformed erp? 1=yes, 0=no
CSD=1;
close all
clc
addpath(genpath('functions/'));
addpath(genpath('../CSDtoolbox/'));
addpath(genpath('../eeglab13_6_5b/'));
%eeglab
chanlocs = readlocs ('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP with reference channel 'FCz' included - hence 65 chans

%% Create paths for subjects
%%
subject_folder = {'AB99C','AG26C','AH05D','AL20D','AM23D','AZ60C','BB12D','BB47C','BD38D','BK48C',...
    'BL31D','CA30D','CC27D','CC76C','CD11D','CF80C','CG52C','CH78C','DN64C','DR69C',...
    'DR79C','DW45D','EB68D','ED62C','EF65C','EH43D','EV61C','FH003D','FO37D','GK67C',...
    'HT17D','HT72C','IA14D','IM19D','JC002D','JD35C','JD40C','JE04D','JE18D','JG10D',...
    'JG77C','JH63C','JK33D','JM25C','JM53D','JN21C','JP82C','JS24D','JS73C','JW71C',...
    'LB49C','LK32D','MB13D','MB16D','MD36C','MG01D','MK84C','MW70C','MY83C','NP44C',...
    'NV75C','OB50C','OH08D','OK34D','PD06D','RD15D','SH81C','SS85C','TB28C','TG51C',...
    'TK66C','VR86C','WH09D','WM22D','ZB41D'};
allsubj =       {'AB99CT','AG26CT','AH05DT','AL20DT','AM23DT','AZ60CT','BB12DT','BB47CT','BD38DT','BK48CT',...
    'BL31DT','CA30DT','CC27DT','CC76CT','CD11DT','CF80CT','CG52CT','CH78CT','DN64CT','DR69CT',...
    'DR79CT','DW45DT','EB68DT','ED62CT','EF65CT','EH43DT','EV61CT','FH003DT','FO37DT','GK67CT',...
    'HT17DT','HT72CT','IA14DT','IM19DT','JC02DT','JD35CT','JD40CT','JE04DT','JE18DT','JG10DT',...
    'JG77CT','JH63CT','JK33DT','JM25CT','JM53DT','JN21CT','JP82CT','JS24DT','JS73CT','JW71CT',...
    'LB49CT','LK32DT','MB13DT','MB16DT','MD36CT','MG01DT','MK84CT','MW70CT','MY83CT','NP44CT',...
    'NV75CT','OB50CT','OH08T','OK34DT','PD06DT','RD15DT','SH81CT','SS85CT','TB28CT','TG51CT',...
    'TK66CT','VR86CT','WH09DT','WM22DT','ZB41DT'};

Dyslexic = {'AH05D','AL20D','AM23D','BB12D','BD38D','BL31D','CA30D','CC27D','CD11D','DW45D'...
    'EB68D','EH43D','FH003D','FO37D','HT17D','IA14D','IM19D','JC002D','JE04D','JE18D'...
    'JG10D','JK33D','JM53D','JS24D','LK32D','MB13D','MB16D','MG01D','OH08D','OK34D'...
    'PD06D','RD15D','WH09D','WM22D','ZB41D'};
AM_control = {'AB99C','AG26C','BK48C','CF80C','DN64C','DR79C','ED62C','EV61C','GK67C','HT72C'...
    'JD35C','JD40C','JG77C','JH63C','JM25C','JN21C','JP82C','LB49C','MW70C','NP44C'...
    'OB50C','SH81C','SS85C'};
RM_control = {'AZ60C','BB47C','CC76C','CG52C','CH78C','DR69C','EF65C','JS73C','JW71C','MD36C'...
    'MK84C','MY83C','NV75C','TB28C','TG51C','TK66C','VR86C'};
%% Define TCD and Monash subjects
TCD_bigdots = {};
Monash_bigdots =  {'AB99CT','AG26CT','AH05DT','AL20DT''AM23DT','AZ60CT','BB12DT','BB47CT','BD38DT','BK48CT',...
    'BL31DT','CA30DT','CC27DT','CC76CT','CD11DT','CF80CT','CG52CT','CH78CT','DN64CT','DR69CT',...
    'DR79CT','DW45DT','EB68DT','ED62CT','EF65CT','EH43DT','EV61CT','FH003DT','FO37DT','GK67CT',...
    'HT17DT','HT72CT','IA14DT','IM19DT','JC02DT','JD35CT','JD40CT','JE04DT','JE18DT','JG10DT',...
    'JG77CT','JH63CT','JK33DT','JM25CT','JM53DT','JN21CT','JP82CT','JS24DT','JS73CT','JW21CT',...
    'LB49CT','LK32DT','MB13DT','MB16DT','MD36CT','MG01DT','MK84CT','MW70CT','MY83CT','NP44CT',...
    'NV75CT','OB50CT','OH08DT','OK34DT','PD06DT','RD15DT','SH81CT','SS85CT','TB28CT','TG51CT',...
    'TK66CT','VR86CT','WH09DT','WM22DT','ZB41DT'};

%%
% DAT1_tags = {'0/1 DAT1 10-repeats','2 DAT1 10-repeats'};
side_tags = {'SLHL','SLHR','SRHL','SRHR'};

%%
duds = []; %% 1(LK_07_04_14) completed the wrong paradigm
single_participants = [subj];
%%
if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
end


%% Define channels, having combined Brain Products and Biosemi data
plot_chans = 1:65;
% exclude_chans = [];
%plot_chans = [1:64];
left_hemi = [1,33,34,3,37,4,38,9,43,8,42,41,12,47,13,48,19,52,18,51,17,23,56,24,57,61,60,28,29];
right_hemi = [2,35,36,39,6,40,7,10,44,11,45,46,49,15,50,16,20,54,21,55,22,58,26,59,27,63,64,31,32];
centre_chans = [5,65,14,53,25,62,30];
elec_pairs = [1,2;33,36;34,35;3,7;37,40;4,6;38,39;...
    41,46;42,45;8,11;43,44;9,10;48,49;13,15;47,50;12,16;
    17,22;51,55;18,21;52,54;19,20;...
    23,27;56,59;24,26;57,58;60,64;61,63;29,31;28,32];

tester = zeros(65,1);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','numbers','plotchans',1:65);
figure
topoplot(tester,chanlocs,'maplimits', ...
    [min(tester)  max(tester)],'electrodes','labels','plotchans',1:65);
%% Triggers

% ITI,left/right
targcodes = zeros(3,2);
targcodes(1,:) = [101,102];
targcodes(2,:) = [103,104];
targcodes(3,:) = [105,106];

%%
fs=500;
numch=65;
rtlim=[0.300 2.100];

ch_contra = [right_hemi;left_hemi];
ch_ipsi = [left_hemi;right_hemi];

allchans = [elec_pairs(:,1)', centre_chans,elec_pairs(:,2)'];
hemi{1} = [elec_pairs(:,2)', centre_chans, elec_pairs(:,1)']; % all contra , center , all ipsi, left patch/lefthand
hemi{2} = [elec_pairs(:,1)', centre_chans, elec_pairs(:,2)']; % all contra , center , all ipsi, right patch/righthand
%%
% stim-locked erps64
ts =  -1*fs:2*fs;
t = ts*1000/fs;

ts_crop = -0.500*fs:1.500*fs;
t_crop = ts_crop*1000/fs;


% pupil uses another time
ts_pupil = -1*fs:4*fs;
t_pupil=ts_pupil*1000/fs;
trs_pupil = [-0.7*fs:fs*2];
tr_pupil = trs_pupil*1000/fs;

% resp-locked erps
trs = [-.700*fs:fs*0.100];
tr = trs*1000/fs;

BL_erp = [-100,0];
BL_pupil=[-100,0];
BL_beta = [-100 0];
BL_alpha = [-100 0];


% zscore threshold
z_thresh = 3;
% define the channels
ch_betai{1} =[43 13]; ch_betai{2} = [15];
ch_betac{1} =[15]; ch_betac{2} = [13 43];
ch_lr{1}= [27]; ch_lr{2}= [23];
ch_rl{1} = [23]; ch_rl{2} = [27];
ch_CPP=[25 53];
ch_alpha{1} = [61]; ch_alpha{2} =[63] ;
ch_alphaC{1}=ch_alpha{2};ch_alphaC{2}=ch_alpha{1};
ch_alphaI{1}=ch_alpha{1};ch_alphaI{2}=ch_alpha{2};

%% Start loop
nameUsed='HAPPE';
fileUsed = [num2str(LPFused) 'fft'];
path_temp = ['../Data/ProcessedSubjects/' nameUsed ]; %test on external HD
for s=1:length(allsubj)
    
    load([path_temp '/' dyscAMRM '/' allsubj{s} '_' fileUsed  'alldots.mat'],'erp','erp_CSD','erp_HPF','allRT','allrespLR','allTrig',...
        'artifact_BLintTo100msPostResponse_n','hand_used','rejected_trial_n','falsealarm','Pupil','ET_BL_resp_artrej');
    % beta times
    STFT_time=[];
    no_of_cycles = 8;
    fs_STFT = 15:25;%good delta range: 8:33, beta is 15:25
    stftlen_STFT = round((1000/round(median(fs_STFT))*no_of_cycles)/2);
    skip_step_beta = 10;
    cc=1;
    for tt = 1:skip_step_beta:length(ts)-(stftlen_STFT)
        tf = tt:tt+stftlen_STFT-1;
        nfft = length(tf);
        freq_temp_STFT = (0:ceil((nfft+1)/2)-1)*fs/nfft;
        STFT_time(cc) = mean(t(tf));
        cc=cc+1;
    end
    
    pause(1)
    
    
    % create Beta waveforms
    if CSD
        erp_beta=double(erp_CSD); %35Hz since looking at beta
        erpalpha = erp_beta;
        erp=erp_CSD;
    else
        erp_beta=double(erp); %35Hz since looking at beta
        erpalpha = erp_beta;
        erp=erp;
    end
    
    baseline_erp = mean(erp_beta(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp_beta = erp_beta-repmat(baseline_erp,[1,size(erp_beta,2),1]); % baseline full erp
    disp('Calculating STFT...')
    STFT = [];
    for trial = 1:size(erp_beta,3)
        cc=1;
        for tt = 1:skip_step_beta:size(erp_beta,2)-(stftlen_STFT)
            tf = tt:tt+stftlen_STFT-1; % define time window
            ep = squeeze(erp_beta(:,tf,trial)); % chop out chan x time window
            nfft = size(ep,2);
            ep = detrend(ep')'; % detrend
            fftx = abs(fft(ep,[],2))./(stftlen_STFT/2);
            fftx = fftx(:,1:ceil((nfft+1)/2));
            ind = find(freq_temp_STFT>fs_STFT(1) & freq_temp_STFT<fs_STFT(end) & freq_temp_STFT~=21);%exclude 21Hz since SSVEP
            %[~,ind] = min(abs(freq_temp_STFT-fs_STFT)); % if you want SSVEP at one particular frequency
            STFT(:,cc,trial) = mean(fftx(:,ind),2);
            cc=cc+1;
        end
        [~,targetSide] = find(allTrig(trial)==targcodes);
        % asymmetry: right minus left. more positive = more right hemi alph
        beta_asym(right_hemi,:,trial) = (STFT(ch_contra(targetSide,:),:,trial)-STFT(ch_ipsi(targetSide,:),:,trial))./...
            ((STFT(ch_contra(targetSide,:),:,trial)+STFT(ch_ipsi(targetSide,:),:,trial))/2);
    end
    
    % create alpha waveforms
    disp('Calculating pre-target Alpha...')
    alpha_bandlimits = [6 11]; % defining the filter for alpha bandpass.
    [H,G]=butter(4,[2*(alpha_bandlimits(1)/fs) 2*(alpha_bandlimits(2)/fs)]); % alpha bandpass for 500Hz

    window = 20; % in samples. Time is double this.
    skip_step_alpha = window/2;

    % Alpha time
    alpha_t=[]; cca=1;
    for tt = 1:skip_step_alpha:length(t_crop)-window
        alpha_t(:,cca) = mean(t_crop(tt:tt+window-1));
        cca=cca+1;
    end
    % Alpha Spectrotemporal Evolution a la Thut
    alpha_TSE = []; alpha_asym = [];
    for trial = 1:size(erp,3)
        % filtering to alpha
        ep_filt = filtfilt(H,G,squeeze(erp(:,:,trial))')';
        % chop off ends and rectify
        ep_filt = abs(ep_filt(:,find(t>=t_crop(1) & t<=t_crop(end))));
        % smooth
        cca=1;
        for tt = 1:skip_step_alpha:size(ep_filt,2)-window
            alpha_TSE(:,cca,trial) = mean(ep_filt(:,tt:tt+window-1),2);
            cca=cca+1;
        end
        alpha_asym(:,:,trial)=0*alpha_TSE(:,:,trial);
        [~,targetSide] = find(allTrig(trial)==targcodes);
        % asymmetry: right minus left. more positive = more right hemi alph
        alpha_asym(right_hemi,:,trial) = (alpha_TSE(ch_contra(targetSide,:),:,trial)-alpha_TSE(ch_ipsi(targetSide,:),:,trial))./...
            ((alpha_TSE(ch_contra(targetSide,:),:,trial)+alpha_TSE(ch_ipsi(targetSide,:),:,trial))/2);
    end
    %  
    % Baseline alpha
    baseline_alpha = mean(alpha_TSE(:,find(alpha_t>BL_alpha(1) & alpha_t<=BL_alpha(2)),:),2);
    alpha_TSE_base = alpha_TSE-repmat(baseline_alpha,[1,size(alpha_TSE,2),1]); % baseline full erp
   
    % Baseline beta
    baseline_beta = mean(STFT(:,STFT_time>BL_beta(1) & STFT_time<=BL_beta(2),:),2);
    beta_TSE_base = STFT-repmat(baseline_beta,[1,size(STFT,2),1]); % baseline full erp
    
    % Baseline erp
    baseline_erp = mean(erp(:,find(t>=BL_erp(1) & t<=BL_erp(2)),:),2);
    erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp
    
    disp(['Subject ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(length(find(allTrig)))])
    
    % Make reponselocked waveforms
    % for each trial, shift the erp signature to make RT = 0;
    alpha_tr= -400:skip_step_alpha*2:100;
    %Response locked STFT time in samples
    alpha_trs = -0.4/(skip_step_alpha/fs):.100/(skip_step_alpha/fs);
    STFT_timer= -400:skip_step_beta*2:100;
    %Response locked STFT time in samples
    STFT_timers = -0.4/(skip_step_beta/fs):.100/(skip_step_beta/fs);
    
    erpr = zeros(size(erp,1),length(tr),size(erp,3));
    Pupilr = zeros(size(Pupil,1),length(tr_pupil),size(Pupil,3));
    STFTr = zeros(size(erp,1),length(STFT_timer),size(erp,3));
    alphar = zeros(size(erp,1),length(alpha_tr),size(erp,3));
    alphar_asym= zeros(size(erp,1),length(alpha_tr),size(erp,3));
    
    validrlock = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT)
        [blah,RTsamp] = min(abs(t*fs/1000-allRT(n)));
        [blah,RTsampBeta] = min(abs(STFT_time*fs/1000-allRT(n))); % get the sample point of the RT.
        [blah,RTsampAlpha] = min(abs(alpha_t*fs/1000-allRT(n)));
        if      RTsamp+trs(1) >0 & RTsamp+trs(end)<=length(t) &...
                RTsampBeta+STFT_timers(1) >0 & RTsampBeta+STFT_timers(end)<=length(STFT_time) & allRT(n)>0 & ...
                RTsampAlpha+alpha_trs(1) >0 & RTsampAlpha+alpha_trs(end)<=length(alpha_t)  % is the RT larger than 1st stim RT point, smaller than last RT point.
            
            erpr(:,:,n) = erp(:,RTsamp+trs,n);
            STFTr(:,:,n) = beta_TSE_base(:,RTsampBeta+STFT_timers,n);
            alphar(:,:,n) = alpha_TSE_base(:,RTsampAlpha+alpha_trs,n);
            alphar_asym(:,:,n) = alpha_asym(:,RTsampAlpha+alpha_trs,n);
            
            validrlock(n)=1;
        end
    end
    
    % create the trials
    for side=1:2
        for hand = 1:2
            for iti = 1:3
                conds{s,iti,side,hand} = find(allTrig==targcodes(iti,side) & allrespLR==1 & ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs & ~artifact_BLintTo100msPostResponse_n & validrlock & ~rejected_trial_n  &...
                    hand_used==hand);
                condsfalse{s,iti,side,hand} = find(allTrig==targcodes(iti,side) & hand_used == hand  & ~artifact_BLintTo100msPostResponse_n);
                RTs{s,iti,side,hand} = allRT([conds{s,iti,side,hand}])*1000/fs;
                RTs_log{s,iti,side,hand} = log(allRT([conds{s,iti,side,hand}])*1000/fs);
                RT_zs{s,iti,side,hand} = zscore([RTs_log{s,iti,side,hand}]);
                RT_factors(s,iti,side,hand) = mean([RTs{s,iti,side,hand}]); % can't include AR_08_04_14 & MH_14_04_14 because of mistake
                
                hit{s,iti,side,hand} = find(allTrig==targcodes(iti,side) & allrespLR==1 & ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs& hand_used==hand);
                miss{s,iti,side,hand} = find(allTrig==targcodes(iti,side) & allrespLR==3 & hand_used==hand);
                wrong{s,iti,side,hand} = find(allTrig==targcodes(iti,side) & allrespLR==2 & hand_used==hand);
                allTs{s,iti,side,hand} = find(allTrig==targcodes(iti,side) &  ...
                    allRT>rtlim(1)*fs & allRT<rtlim(2)*fs& hand_used==hand);
                fAALL{s,iti,side,hand} = falsealarm([condsfalse{s,iti,side,hand}]);
            end
            RT_all(s,side,hand) = mean([RTs{s,:,side,hand}]);
            RT_var(s) = var([RTs{s,:}]);
            fA_all(s,side,hand) = sum([fAALL{s,:,side,hand}]);
            fAtemp = [fAALL{s,:,side,hand}];fAtemp(fAtemp>1)=1;
            fA_count(s,side,hand) = sum(fAtemp);
            temp = [RTs{s,:,side,hand}]; tempz = [RT_zs{s,:,side,hand}];
            temp = temp(find(tempz<z_thresh));
            RT_all_zs(s,side,hand) = mean(temp);
            RT_median_all(s,side,hand) = median([RTs{s,:,side,hand}]);
            RT_log_all(s,side,hand) = mean([RTs_log{s,:,side,hand}]);
            
            hit_all(s,side,hand) = length([hit{s,:,side,hand}]);
            miss_all(s,side,hand) = length([miss{s,:,side,hand}]);
            wrong_all(s,side,hand)=length([wrong{s,:,side,hand}]);
            hit_rate(s,side,hand) = 100*hit_all(s,side,hand)/(miss_all(s,side,hand)+hit_all(s,side,hand)+wrong_all(s,side,hand));
        end
    end
    
    disp(['Subject ',allsubj{s},' Total Valid Trials: ',num2str(length([conds{s,:,:}])), ...
        ' = ',num2str(round(100*(length([conds{s,:,:}]))/(length(find(allTrig))))),'%'])
    trialinfo.accepted = length([conds{:}]);
    trialinfo.alltrials=length(allTrig);
    for side = 1:2
        for hand=1:2
            ERP_side(s,:,:,side,hand) = squeeze(mean(erp(1:numch,:,[conds{s,:,side,hand}]),3));
            ERPr_side(s,:,:,side,hand) = squeeze(mean(erpr(1:numch,:,[conds{s,:,side,hand}]),3));
            
            ERP_hemi_side(s,allchans,:,side,hand) = squeeze(mean(erp(hemi{side},:,[conds{s,:,side,hand}]),3));
            ERPr_hemi_side(s,allchans,:,side,hand) = squeeze(mean(erpr(hemi{side},:,[conds{s,:,side,hand}]),3));
            
            CPP_side(s,:,side,hand) = squeeze(mean(mean(erp(ch_CPP,:,[conds{s,:,side,hand}]),1),3));
            CPPr_side(s,:,side,hand) = squeeze(mean(mean(erpr(ch_CPP,:,[conds{s,:,side,hand}]),1),3));
            N2c_side(s,:,side,hand) = squeeze(mean(mean(erp(ch_lr{side},:,[conds{s,:,side,hand}]),1),3));
            N2i_side(s,:,side,hand) = squeeze(mean(mean(erp(ch_rl{side},:,[conds{s,:,side,hand}]),1),3));
            bERPr_side_all(s,:,:,side,hand) = squeeze(mean(STFTr(1:numch,:,[conds{s,:,side,hand}]),3)); %(s, chan, STFT_timers, targetside)
            bERP_side_all(s,:,:,side,hand) = squeeze(mean(beta_TSE_base(1:numch,:,[conds{s,:,side,hand}]),3));
            
            bERP_r_hemi_side(s,allchans,:,side,hand) = squeeze(mean(STFTr(hemi{hand},:,[conds{s,:,side,hand}]),3)); %(s, chan, STFT_timers, targetside)
            bERP_hemi_side(s,allchans,:,side,hand) = squeeze(mean(beta_TSE_base(hemi{hand},:,[conds{s,:,side,hand}]),3));
            
            betai_side(s,:,side,hand) = squeeze(mean(mean(beta_TSE_base(ch_betai{hand},:,[conds{s,:,side,hand}]),1),3));
            betac_side(s,:,side,hand) = squeeze(mean(mean(beta_TSE_base(ch_betac{hand},:,[conds{s,:,side,hand}]),1),3));
            betarc_side(s,:,side,hand) = squeeze(mean(mean(STFTr(ch_betac{hand},:,[conds{s,:,side,hand}]),1),3));
            betari_side(s,:,side,hand) = squeeze(mean(mean(STFTr(ch_betai{hand},:,[conds{s,:,side,hand}]),1),3));
            
            alpha_side(s,allchans,:,side,hand) = squeeze(mean(alpha_TSE(hemi{side},:,[conds{s,:,side,hand}]),3));
            alpha_base_side(s,allchans,:,side,hand) = squeeze(mean(alpha_TSE_base(hemi{side},:,[conds{s,:,side,hand}]),3)); %(s, chan, STFT_timers, targetside)
            alpha_asym_side(s,:,:,side,hand) = squeeze(mean(alpha_asym(1:numch,:,[conds{s,:,side,hand}]),3));
            
            alpha_asym_avg_side(s,right_hemi,:,side) = (alpha_side(s,right_hemi,:,side)-alpha_side(s,left_hemi,:,side))./...
                ((alpha_side(s,right_hemi,:,side)+alpha_side(s,left_hemi,:,side))/2);
            
            alphar_base_side(s,allchans,:,side,hand) = squeeze(mean(alphar(hemi{side},:,[conds{s,:,side,hand}]),3)); %(s, chan, STFT_timers, targetside)
            alphar_asym_side(s,:,:,side,hand) = squeeze(mean(alphar_asym(1:numch,:,[conds{s,:,side,hand}]),3));
            
            %% Code adapted from Ger's Current Biology cpp code to pull out CPP onset latency, also use it for beta onset latency
            % Define CPP onset search window, from 0 to 1000ms
            CPP_search_t  = [0,1000];beta_search_t=[0,1500];
            % Same window in samples
            CPP_search_ts  = [find(t==CPP_search_t(1)),find(t==CPP_search_t(2))];
            % Size of sliding window. This is in fact 1/4 of the search window in ms.
            % So 25 is 100ms. (25 samples x 2ms either side of a particular sample).
            max_search_window = 25;consecutive_windows=10;
            CPP_temp = squeeze(mean(erp(ch_CPP,:,[conds{s,:,side,hand}]),1)); % time x trial
            CPPs(:,side,hand) = squeeze(mean(CPP_temp(:,:),2)); % average across trial for plot later on, not used to find onsets.
            CPP_side_onsets(s,side,hand) = getOnsetfn(CPP_temp,CPP_search_t,CPP_search_ts,max_search_window,consecutive_windows,t,subject_folder(s));
            
            CPPampvar = squeeze(erpr(ch_CPP,find(tr==0),[conds{s,:,:,:}]));
            
            %contral and ipsi betas
            beta_temp =  squeeze(mean(beta_TSE_base(ch_betac{hand},:,[conds{s,:,side,hand}]),1)); % time x trial
            betasc(:,side,hand) = squeeze(mean(beta_temp(:,:),2));
            betac_side_onsets(s,side,hand) = getOnsetfn(beta_temp,beta_search_t,CPP_search_ts,25,2,STFT_time,subject_folder(s));%smaller windows due to large time skips
            
            clear beta_temp
            beta_temp =  squeeze(mean(beta_TSE_base(ch_betai{hand},:,[conds{s,:,side,hand}]),1)); % time x trial
            betasi(:,side,hand) = squeeze(mean(beta_temp(:,:),2));
            betai_side_onsets(s,side,hand) = getOnsetfn(beta_temp,beta_search_t,CPP_search_ts,25,2,STFT_time,subject_folder(s));%smaller windows due to large time skips
            clear beta_temp
        end
    end

    %%   plot topology of ERPr
    close all
    t1 = -300; t2 = 100;
    plot_mean = squeeze(nanmean(mean(mean(ERPr_hemi_side(1,:,find(tr>=t1 & tr<t2),:,:),3),4),5));
    h1=figure(1)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','numbers','plotchans',plot_chans);
    title([subject_folder{s} ' Percentage ', num2str(round(100*(length([conds{s,:,:}]))/(length(find(allTrig))))),'%'])
    saveas(h1,['Figures/TopoPlotCPP/topoplot', subject_folder{s} '.jpg'])
    h2=figure(2)
    ERPr_plot = squeeze(mean(ERPr_hemi_side,5));
    plottopo(ERPr_plot,'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
        min(min(min(ERPr_plot(:))))  max(max(max(ERPr_plot(:))))], ...
        'title',[subject_folder{s} ' Percentage ', num2str(round(100*(length([conds{s,:,:}]))/(length(find(allTrig))))),'%']...
        ,'legend',side_tags,'showleg','on','ydir',1)
    saveas(h2,['Figures/plottopoCPP/plottopo', subject_folder{s} '.jpg'])
    %%   plot topology of beta
    close all
    t1 = -300; t2 = 100;
    plot_mean = squeeze(nanmean(nanmean(nanmean(bERPr_side_all(1,:,find(STFT_timer>=t1 & STFT_timer<t2),:,1),3),4),5));
    h1=figure(1)
    topoplot(plot_mean,chanlocs,'maplimits', ...
        [min(plot_mean) max(plot_mean)], ...
        'electrodes','numbers','plotchans',plot_chans);
    title([subject_folder{s} ' Percentage ', num2str(round(100*(length([conds{s,:,:}]))/(length(find(allTrig))))),'%'])
    saveas(h1,['Figures/TopoPlotbeta/topoplot', subject_folder{s} '.jpg'])
    h2=figure(2)
    betar_plot = squeeze(mean(bERPr_side_all(1,:,:,:,:),4));
    plottopo(betar_plot,'chanlocs',chanlocs,'limits',[STFT_timer(1) STFT_timer(end) ...
        min(min(min(betar_plot(:))))  max(max(max(betar_plot(:))))], ...
        'title',[subject_folder{s} ' Percentage ', num2str(round(100*(length([conds{s,:,:}]))/(length(find(allTrig))))),'%']...
        ,'legend',side_tags,'showleg','on','ydir',1)
    saveas(h2,['Figures/plottopobeta/plottopo', subject_folder{s} '.jpg'])
    %% plot the time-frequency component of all erps
    erphand1 = erp(:,1:1501,[conds{s,:,:,1}]);
    erphand2 = erp(:,1:1501,[conds{s,:,:,2}]);
    erprhand1 = erpr(:,:,[conds{s,:,:,1}]);
    erprhand2 = erpr(:,:,[conds{s,:,:,2}]);
    %%
    red=[1 0 0];blue=[0 0 1];yellow=[1 1 0];green =[0 1 0];
    colour = [red;blue;yellow;green];
    legend_tags = {'SLHL','SLHR','SRHL','SRHR'};
    clear h
    figure
    hold on
    for side=1:2
        for hand=1:2
            c=(side-1)*2+hand;
            h(c) = plot(tr,squeeze(ERPr_hemi_side(1,25,:,side,hand)),'LineWidth',3,'LineStyle','-','Color',colour(c,:));hold on
            set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
            ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
            xlabel('Time (ms)','FontName','Arial','FontSize',16)
            title(['CPP (resp-locked) ' 'Side' num2str(side) 'Hand' num2str(hand)])
            line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
            line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
        end
    end
    legend(h,legend_tags, ...
        'FontSize',16,'Location','NorthWest');
    %%   plot CPP with onset marked
    h1 = figure;
    colors = {'b' 'r'; 'g' 'm'};
    for side = 1:2
        for hand=1:2
            plot(t,squeeze(CPPs(:,side,hand)),'Color',colors{side,hand},'LineWidth',2), hold on
            line([mean(CPP_side_onsets(s,side,hand),1),mean(CPP_side_onsets(s,side,hand),1)],ylim,'Color',colors{side,hand},'LineWidth',1.5);
            line([0,0],ylim,'Color','k','LineWidth',1);
            line(xlim,[0,0],'Color','k','LineWidth',1);
        end
    end
    title(subject_folder{s})
    saveas(h1,['Figures/CPPonset/CPPonset', subject_folder{s} '.jpg'])
    
    %% Extract Response Locked CPP slope and press level
    % press level
    for side=1:2
        for hand=1:2
            CPP_click(s,side,hand) = squeeze(CPPr_side(s,tr==0,side,hand));
        end
    end
    %CPP build-up defined as the slope of a straight line fitted to the
    %response-locked waveform at during "slope_timeframe_index" defined for
    %each participant here:
    clear slope_timeframe_index
    slope_timeframe_index(s,2)=find(nanmean(nanmean(CPPr_side(s,:,:,:),3),4)==max(mean(mean(CPPr_side(s,find(tr<0),:,:),3),4)));%max amplitude index
    slope_timeframe_index(s,3)=find(tr==0);
    if slope_timeframe_index(s,2)>150
        slope_timeframe_index(s,1)=slope_timeframe_index(s,2)-150;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    else
        slope_timeframe_index(s,1)=1;%move forward 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
    %Now find and save CPPr slope
    for side = 1:2
        for hand=1:2
            coef = polyfit(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,3)),(CPPr_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,3),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            CPPr_slope(s,side,hand)=coef(1);
        end
    end
    CPPall=[];
    for side=1:2
        for hand=1:2
            CPPtimetemp = [];ks=1;
            for i = [conds{s,:,side,hand}]
                CPPsingletrial = mean(erpr(ch_CPP,:,i),1);
                CPPtimetemp(ks) =  ...
                    tr(max(find(CPPsingletrial==max(CPPsingletrial(:,find(tr<=0 & tr>-300))))));
                ks=ks+1;
            end
            CPPr_peakTime(s,side,hand)= mean(CPPtimetemp);
            CPPall = [CPPall CPPtimetemp];
        end
    end
    
    CPPr_peakTimeVar = CPPall;
    %%%Plot each individual participant's CPPr_slope with time-window varying
    %%%per participant
    clear h h1
    h1 = figure;
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(tr,CPPr_side(s,:,side,hand),'LineWidth',3,'LineStyle','-');hold on
            coef = polyfit(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(CPPr_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            r = coef(1) .* tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');hold on;
            line([-50 50], [CPP_click(s,side,hand),CPP_click(s,side,hand)], 'Linewidth',5, 'LineStyle', '-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' CPP (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/CPPAll/CPPAll', subject_folder{s} '.jpg'])
    
    %% extract reponse of beta slope and amplitude
    % betas with onset marked
    clear h1
    h1 = figure;
    colors = {'b' 'r'; 'g' 'm'};
    for side = 1:2
        for hand=1:2
            plot(STFT_time,squeeze(betasc(:,side,hand)),'Color',colors{side,hand},'LineWidth',2), hold on
            line([mean(betac_side_onsets(s,side,hand),1),mean(betac_side_onsets(s,side,hand),1)],ylim,'Color',colors{side,hand},'LineWidth',1.5);
            line([0,0],ylim,'Color','k','LineWidth',1);
            line(xlim,[0,0],'Color','k','LineWidth',1);
        end
    end
    title(subject_folder{s})
    saveas(h1,['Figures/BetaOnset/BetaOnset', subject_folder{s} '.jpg'])
    % beta with amplitude and slopes
    for side=1:2
        for hand=1:2
            Betac_click(s,side,hand) = squeeze(betarc_side(s,STFT_timer==0,side,hand));
            Betai_click(s,side,hand) = squeeze(betari_side(s,STFT_timer==0,side,hand));
        end
    end
    %Beta build-up defined as the slope of a straight line fitted to the
    %response-locked waveform at during "slope_timeframe_index" defined for
    %each participant here:
    clear slope_timeframe_index
    slope_timeframe_index(s,2)=find(nanmean(nanmean(-betarc_side(s,:,:,:),3),4)==max(mean(mean(-betarc_side(s,find(STFT_timer>=-100 & STFT_timer<0),:,:),3),4)));%max amplitude index
    if slope_timeframe_index(s,2)>15
        slope_timeframe_index(s,1)=slope_timeframe_index(s,2)-15;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    else
        slope_timeframe_index(s,1)=1;%move forward 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
    %Now find and save CPPr slope
    for side = 1:2
        for hand=1:2
            coef = polyfit(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(betarc_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            betac_slope(s,side,hand)=coef(1);
        end
    end
    
    %betaislope
    clear slope_timeframe_index
    slope_timeframe_index(s,2)=find(mean(mean(-betari_side(s,:,:,:),3),4)==max(mean(mean(-betari_side(s,find(STFT_timer>=-100 & STFT_timer<0),:,:),3),4)));%max amplitude index
    if slope_timeframe_index(s,2)>10
        slope_timeframe_index(s,1)=slope_timeframe_index(s,2)-10;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    else
        slope_timeframe_index(s,1)=1;%move forward 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
    %Now find and save CPPr slope
    for side = 1:2
        for hand=1:2
            coef = polyfit(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(betari_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            betai_slope(s,side,hand)=coef(1);
        end
    end
    %% Plot each individual participant's betar_slope with time-window varying
    %%%per participant
    clear h h1
    h1 = figure;
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(STFT_timer,betari_side(s,:,side,hand),'LineWidth',3,'LineStyle','-');hold on
            coef = polyfit(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(betari_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            r = coef(1) .* STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');hold on;
            line([-50 50], [Betai_click(s,side,hand),Betai_click(s,side,hand)], 'Linewidth',5, 'LineStyle', '-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' Beta (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/BetaAll/BetaiAll', subject_folder{s} '.jpg'])
    
    %%%Plot each individual participant's betar_slope with time-window varying
    %%%per participant
    clear h h1
    h1 = figure;
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(STFT_timer,betarc_side(s,:,side,hand),'LineWidth',3,'LineStyle','-');hold on
            coef = polyfit(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(betarc_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            r = coef(1) .* STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(STFT_timer(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');hold on;
            line([-50 50], [Betac_click(s,side,hand),Betac_click(s,side,hand)], 'Linewidth',5, 'LineStyle', '-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' Beta (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/BetaAll/BetacAll', subject_folder{s} '.jpg'])
    
%% extract reponse of alpha slope and amplitude
    alpha_c_side(s,:,:,:) = squeeze(mean(alphar_base_side(s,ch_alpha{1},:,:,:),2));
    alpha_i_side(s,:,:,:) = squeeze(mean(alphar_base_side(s,ch_alpha{2},:,:,:),2));
    % alpha with amplitude and slopes
    for side=1:2
        for hand=1:2
            alphac_click(s,side,hand) = squeeze(alpha_c_side(s,alpha_tr==0,side,hand));
            alphai_click(s,side,hand) = squeeze(alpha_i_side(s,alpha_tr==0,side,hand));
        end
    end    
    %alpha build-up defined as the slope of a straight line fitted to the
    %response-locked waveform at during "slope_timeframe_index" defined for
    %each participant here:
    clear slope_timeframe_index
    slope_timeframe_index(s,2)=find(mean(mean(-alpha_c_side(s,:,:,:),3),4)==max(mean(mean(-alpha_c_side(s,find(alpha_tr==0),:,:),3),4)));%max amplitude index
    if slope_timeframe_index(s,2)>8
        slope_timeframe_index(s,1)=slope_timeframe_index(s,2)-8;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    else
        slope_timeframe_index(s,1)=1;%move forward 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
    %Now find and save alphac slope
    for side = 1:2
        for hand=1:2
            coef = polyfit(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(alpha_c_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            alphac_slope(s,side,hand)=coef(1);
        end
    end
    
    %alphaislope
    clear slope_timeframe_index
    slope_timeframe_index(s,2)=find(mean(mean(-alpha_i_side(s,:,:,:),3),4)==max(mean(mean(-alpha_i_side(s,find(alpha_tr==0),:,:),3),4)));%max amplitude index
    if slope_timeframe_index(s,2)>8
        slope_timeframe_index(s,1)=slope_timeframe_index(s,2)-8;%subtract 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    else
        slope_timeframe_index(s,1)=1;%move forward 50samples (i.e. 100ms) from max amplitude index to form slope_timeframe window
    end
    %Now find and save alphai slope
    for side = 1:2
        for hand=1:2
            coef = polyfit(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(alpha_i_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            alphai_slope(s,side,hand)=coef(1);
        end
    end
    %% Plot each individual participant's alphar_slope with time-window varying
    %%%per participant
    clear h h1
    h1 = figure;
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(alpha_tr,alpha_i_side(s,:,side,hand),'LineWidth',3,'LineStyle','-');hold on
            coef = polyfit(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(alpha_i_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            r = coef(1) .* alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');hold on;
            line([-50 50], [alphai_click(s,side,hand),alphai_click(s,side,hand)], 'Linewidth',5, 'LineStyle', '-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' Alpha Ipsi (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/BetaAll/AlphaiAll', subject_folder{s} '.jpg'])
    
    %%%Plot each individual participant's betar_slope with time-window varying
    %%%per participant
    clear h h1
    h1 = figure;
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(alpha_tr,alpha_c_side(s,:,side,hand),'LineWidth',3,'LineStyle','-');hold on
            coef = polyfit(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)),(alpha_c_side(s,slope_timeframe_index(s,1):slope_timeframe_index(s,2),side,hand)),1);% coef gives 2 coefficients fitting r = slope * x + intercept
            r = coef(1) .* alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)) + coef(2); %r=slope(x)+intercept, r is a vectore representing the linear curve fitted to the erpr during slope_timeframe
            plot(alpha_tr(slope_timeframe_index(s,1):slope_timeframe_index(s,2)), r,'Linewidth',2, 'LineStyle', ':');hold on;
            line([-50 50], [alphac_click(s,side,hand),alphac_click(s,side,hand)], 'Linewidth',5, 'LineStyle', '-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-500,100]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title([subject_folder{s}, ' Alpha Contra (resp-locked) by Hemifield'])
    line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags,'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/BetaAll/AlphacAll', subject_folder{s} '.jpg'])  
    %% Extract N2c and N2i latency and amplitude :
    for side=1:2
        for hand=1:2
            %         to use each participant's average N2c/i to get their peak latency index:
            N2c=N2c_side(s, :, side,hand);
            N2i=N2i_side(s, :, side,hand);
            N2c_peak_amp_index=find(N2c==min(N2c(find(t==150):find(t==450))));%Find Left target max peak latency for N2c
            N2i_peak_amp_index=find(N2i==min(N2i(find(t==150):find(t==450))));%Find Left target max peak latency for N2i
            N2c_peak_amp_index_t(s,side,hand)=t(N2c_peak_amp_index);%Find max peak latency for N2c in ms
            N2i_peak_amp_index_t(s,side,hand)=t(N2i_peak_amp_index);%Find max peak latency for N2i in ms
            %N2 Latency:
            N2cN2i_latency_ByTargetSide(s,side,hand,:) = [N2c_peak_amp_index_t(s,side,hand),N2i_peak_amp_index_t(s,side,hand)]; %(LeftTargetN2c_latency, RightTargetN2c_latency, LeftTargetN2i_latency, RightTargetN2i_latency)
            
            window=25; %this is the time (in samples) each side of the peak latency (so it's 50ms each side of peak latency - so a 100ms window)
            N2c = squeeze(mean(mean(mean(N2c_side,1),3),4)); % time
            N2i = squeeze(mean(mean(mean(N2i_side,1),3),4)); % time
            
            
            max_peak_N2c(s,side,hand)=squeeze(mean(N2c_side(s,N2c_peak_amp_index-window:N2c_peak_amp_index+window, side,hand),2));
            max_peak_N2i(s,side,hand)=squeeze(mean(N2i_side(s,N2i_peak_amp_index-window:N2i_peak_amp_index+window, side,hand),2));
        end
    end
    N2cN2i_amp_ByTargetSide_ParticipantLevel = [max_peak_N2c,max_peak_N2i]; %(LeftTargetN2c, RightTargetN2c, LeftTargetN2i, RightTargetN2i)
    
    %Plot N2c
    clear h
    h1 = figure
    hold on
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(t,squeeze(N2c_side(s,:,side,hand)),'LineWidth',3,'LineStyle','-', 'Color',colors{side,hand});hold on
            line([t(N2c_peak_amp_index-window),t(N2c_peak_amp_index+window)],[max_peak_N2c(s,side,hand), max_peak_N2c(s,side,hand)],'Color',colors{side,hand}, 'Linewidth',2,'LineStyle','-');
        end
    end
    
    set(gca,'FontSize',16,'xlim',[-100,600],'xtick',[-100,0:200:600]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['Subj: ',num2str(s),' N2c by Hemifield'])
    line([0,0],[0 -30],'Color','k','LineWidth',1.5,'LineStyle','--');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags, 'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/N2c/N2c', subject_folder{s} '.jpg'])
    
    %Plot N2i
    clear h
    h1 = figure
    hold on
    for side = 1:2
        for hand=1:2
            h((side-1)*2+hand) = plot(t,squeeze(N2i_side(s,:,side,hand)),'LineWidth',3,'LineStyle','-', 'Color',colors{side,hand});hold on
            line([t(N2i_peak_amp_index-window),t(N2i_peak_amp_index+window)],[max_peak_N2i(s,side,hand), max_peak_N2i(s,side,hand)],'Color',colors{side,hand}, 'Linewidth',2,'LineStyle','-');
        end
    end
    set(gca,'FontSize',16,'xlim',[-100,600],'xtick',[-100,0:200:600]);%,'ylim',[-4,8],'ytick',[-4:2:8]);%,'ylim',[-1.5,0.5]);
    ylabel('Amplitude (\muVolts)','FontName','Arial','FontSize',16)
    xlabel('Time (ms)','FontName','Arial','FontSize',16)
    title(['Subj: ',num2str(s),' N2i by Hemifield'])
    line([0,0],[0 -30],'Color','k','LineWidth',1.5,'LineStyle','-');
    line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
    legend(h,side_tags, 'FontSize',16,'Location','NorthWest');
    pause(1)
    saveas(h1,['Figures/N2i/N2i', subject_folder{s} '.jpg'])
    
    %% Extract Alpha Power (Pre and Post)
    for side = 1:2
        for hand=1:2
            preAlpha(:,side,hand)= squeeze(mean(mean(mean(alpha_TSE([61 62 63],find(alpha_t>-500 & alpha_t<-50),[conds{s,:,side,hand}]),1),2),3));
            preAlpha_L(:,side,hand)=squeeze(mean(mean(mean(alpha_TSE(ch_alphaC{side},find(alpha_t>-500 & alpha_t<-50),[conds{s,:,side,hand}]),1),2),3));
            preAlpha_R(:,side,hand)=squeeze(mean(mean(mean(alpha_TSE(ch_alphaI{side},find(alpha_t>-500 & alpha_t<-50),[conds{s,:,side,hand}]),1),2),3));
            preAlpha_asym(:,side,hand)=squeeze(mean(mean(mean(alpha_asym(ch_alpha{2},find(alpha_t>-500 & alpha_t<-50),[conds{s,:,side,hand}]),1),2),3));
            
            postAlpha(:,side,hand)=squeeze(mean(mean(mean(alpha_TSE_base([61 62 63],find(alpha_t>150 & alpha_t<1000),[conds{s,:,side,hand}]),1),2),3));
            postAlpha_L(:,side,hand)=squeeze(mean(mean(mean(alpha_TSE_base(ch_alphaC{side},find(alpha_t>150 & alpha_t<1000),[conds{s,:,side,hand}]),1),2),3));
            postAlpha_R(:,side,hand)=squeeze(mean(mean(mean(alpha_TSE_base(ch_alphaI{side},find(alpha_t>150 & alpha_t<1000),[conds{s,:,side,hand}]),1),2),3));
            postAlpha_asym(:,side,hand)=squeeze(mean(mean(mean(alpha_asym(ch_alpha{2},find(alpha_t>150 & alpha_t<1000),[conds{s,:,side,hand}]),1),2),3));
            
            postAlphaR(:,side,hand)=squeeze(mean(mean(mean(alphar([61 62 63],find(alpha_tr==0),[conds{s,:,side,hand}]),1),2),3));
            postAlphaR_L(:,side,hand)=squeeze(mean(mean(mean(alphar(ch_alphaC{side},find(alpha_tr==0),[conds{s,:,side,hand}]),1),2),3));
            postAlphaR_R(:,side,hand)=squeeze(mean(mean(mean(alphar(ch_alphaI{side},find(alpha_tr==0),[conds{s,:,side,hand}]),1),2),3));
            postAlphaR_asym(:,side,hand)=squeeze(mean(mean(mean(alphar_asym(ch_alpha{2},find(alpha_tr==0),[conds{s,:,side,hand}]),1),2),3));
            
        end
    end
    %% extract pupil diameter
    pupiltype=1;
    for ks = [2 3 4]
        Pupiltemp = squeeze(Pupil(ks,:,:));%6,
        Pupilrtemp = squeeze(Pupilr(ks,:,:));%8
        
        
        Pupiltemp = diag(1/max(max(Pupiltemp,[],2)))*Pupiltemp;
        Pupilrtemp = diag(1/max(max(Pupiltemp,[],2)))*Pupilrtemp;
        
        baseline  = mean(Pupiltemp(find(t_pupil>=BL_pupil(1) & t_pupil<=BL_pupil(2)),:),1);
        Pupiltemp     = Pupiltemp - repmat(baseline,[size(Pupiltemp,1)],1); % baseline full erp
        Pupilrtemp    = Pupilrtemp - repmat(baseline,[size(Pupilrtemp,1)],1); % baseline full erp
        
        diameteratTarget = mean(Pupiltemp(find(tr_pupil>=-100 & tr_pupil<=100),:),1);
        
        PupilN(pupiltype,:,:) = Pupiltemp;
        PupilrN(pupiltype,:,:) = Pupilrtemp;
        Pupilbaseline(pupiltype,:,:) = baseline;
        Pupildiametertarget(pupiltype,:,:) = diameteratTarget;
        pupiltype = pupiltype+1;
    end
    %% pupil diamteter
    for side=1:2
        for hand=1:2
            PupilD(:,:,side,hand) = mean(PupilN(:,:,[conds{s,:,side,hand}]),3);
            PupilrD(:,:,side,hand)= mean(PupilrN(:,:,[conds{s,:,side,hand}]),3);
        end
    end    
end



%pupil.baselinephase.bp  = squeeze(angle(mean(exp(1i*(inp_pupil(8,find(tt>=set.BL(1) & tt<=set.BL(2)),:)))))); %phase pupil

%% save the data for faster plots
if CSD
    save(['Data/ERPs/' nameUsed '/CSD_filt' num2str(LPFused) '_' dyscAMRM '_group_plots_' num2str(single_participants)],'CPP_side','CPPr_side','chanlocs','ERP_side','ERPr_side',...
        'N2c_side','N2i_side','RTs','t','plot_chans','side_tags','tr','ERP_hemi_side','ERPr_hemi_side','hit','wrong','miss','allTs',...
        'CPP_side_onsets','CPPr_slope','CPP_click','CPPr_peakTime','N2i_peak_amp_index_t','N2c_peak_amp_index_t','max_peak_N2i','max_peak_N2c','RT_all','hit_rate',...
        'betac_side_onsets','betai_side_onsets','betac_slope','betai_slope','Betac_click','Betai_click',...
        'betac_side','betai_side','betarc_side','betari_side','bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side',...
        'STFT_timer','STFT_time','alpha_t','fAALL','fA_all','fA_count',...
        'alpha_side','alpha_base_side','alpha_asym_side','alpha_asym_avg_side','alphar_base_side','alpha_tr','alphar_asym_side',...
        'alphac_click','alphai_click',...
        't_pupil','tr_pupil','PupilD','PupilrD','Pupilbaseline','Pupildiametertarget',...
        'preAlpha','preAlpha_L','preAlpha_R','preAlpha_asym','postAlpha','postAlpha_L','postAlpha_R','postAlpha_asym',...
        'postAlphaR','postAlphaR_L','postAlphaR_R','postAlphaR_asym','alphai_slope','alphac_slope',...
        'RT_var','CPPampvar','CPPr_peakTimeVar','trialinfo','jobID')
    save(['Data/beta_search/' dyscAMRM 'groupdata_CSD',num2str(single_participants)], 'erphand1','erphand2')
else
    save(['Data/ERPs/' nameUsed '/CSD_filt' num2str(LPFused) '_' dyscAMRM '_group_plots_' num2str(single_participants)],'CPP_side','CPPr_side','chanlocs','ERP_side','ERPr_side',...
        'N2c_side','N2i_side','RTs','t','plot_chans','side_tags','tr','ERP_hemi_side','ERPr_hemi_side','hit','wrong','miss','allTs',...
        'CPP_side_onsets','CPPr_slope','CPP_click','CPPr_peakTime','N2i_peak_amp_index_t','N2c_peak_amp_index_t','max_peak_N2i','max_peak_N2c','RT_all','hit_rate',...
        'betac_side_onsets','betai_side_onsets','betac_slope','betai_slope','Betac_click','Betai_click',...
        'betac_side','betai_side','betarc_side','betari_side','bERPr_side_all','bERP_side_all','bERP_r_hemi_side','bERP_hemi_side',...
        'STFT_timer','STFT_time','alpha_t','fAALL','fA_all','fA_count',...
        'alpha_side','alpha_base_side','alpha_asym_side','alpha_asym_avg_side','alphar_base_side','alpha_tr','alphar_asym_side',...
        'alphac_click','alphai_click',...
        't_pupil','tr_pupil','PupilD','PupilrD','Pupilbaseline','Pupildiametertarget',...
        'preAlpha','preAlpha_L','preAlpha_R','preAlpha_asym','postAlpha','postAlpha_L','postAlpha_R','postAlpha_asym',...
        'postAlphaR','postAlphaR_L','postAlphaR_R','postAlphaR_asym','alphai_slope','alphac_slope',...
        'RT_var','CPPampvar','CPPr_peakTimeVar','trialinfo','jobID')
    save(['Data/beta_search/' dyscAMRM 'groupdata_',num2str(single_participants)], 'erphand1','erphand2')
end
return
