targcodes = [101:106];

old_fs = 500; % old sample rate
fs = 500; %new sample rate

%left/right
targcodesLR = zeros(3,2);
targcodesLR(:,1) = [101,103,105];
targcodesLR(:,2) = [102,104,106];

% in sample points, the ERP epoch
ts = -1.0*fs:2*fs; % -1000ms to 1880ms with 200ms either side.
t = ts*1000/fs;
ts_crop = -1*fs:2*fs;
t_crop = ts_crop*1000/fs;


% pupil uses another time
ts_pupil = -1.2*fs:4*fs;
t_pupil=ts_pupil*1000/fs;
ts_pupilcrop = -1*fs:4*fs;
t_pupilcrop = ts_pupilcrop*1000/fs;

BLint = [-100 0];   % baseline interval in ms
default_response_time = 1.1-0.1;

Alpha_samps = ((abs(t_crop(1))+abs(t_crop(end)))/50)-1; %50ms alpha samples
ERP_samps = length(t_crop);

nchan = 65;

LPFcutoff=LPFused;       % Low Pass Filter cutoff
if LPFused>0
    LPF = 1;    % 1 = low-pass filter the data, 0=don't.
else
    LPF=0;
end
if HPFused>0 && HPFused<900
    HPF=1;
    HPFcutoff=HPFused;       % High Pass Filter - either 0.01, 0.1, or 0.25 for cuttoff as I have filters .mats designed for these using "fdatool" tool in MATLAB
else
    HPF=0; %no need for cutoff if using 1Hz original HAPPE protocol or 0 HZ
    HPFcutoff=0;
end



if HPF==1
    if HPFcutoff==0.01;
        load('butter_HPF_0p01');
        HD_HPF = Hd;
    elseif  HPFcutoff==0.1;
        load('butter_HPF_0p1');
        HD_HPF = Hd;
    elseif  HPFcutoff==0.25;
        load('butter_HPF_0p25');
        HD_HPF = Hd;
        
    else
        
        N    = 2;    % Order
        F3dB = HPFcutoff;    % 3-dB Frequency
        Fs   = 500;  % Sampling Frequency
        
        h = fdesign.highpass('n,f3db', N, F3dB, Fs);
        
        HD_HPF = design(h, 'butter');
    end
end

bandlimits(1,1) = 8; % defining the filter for alpha bandpass.
bandlimits(1,2) = 13;

[H1,G1]=butter(4,[2*(bandlimits(1,1)/old_fs) 2*(bandlimits(1,2)/old_fs)]); % alpha bandpass
[H2,G2]=butter(4,[2*(bandlimits(1,1)/fs) 2*(bandlimits(1,2)/fs)]); % alpha bandpass

PretargetARwindow=[-0.500,0];%time window (in seconds, must be factor of the 50ms alpha samples) to search for pre-target artifacts

% frontal channels, occipital channels, and POz and Pz, CPz, CP1, CP2, P1, P2

ARchans = [1:65];
ARchans_for_blinks = [1:65];
artifth = 100;
artifchans=[];  % keep track of channels on which the threshold is exceeded, causing trial rejection

% chanlocs = readlocs('cap128.loc');
chanlocs = readlocs ('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP with reference channel 'FCz' included - hence 65 chans
chanlocs = chanlocs(1:nchan)';

rejected_trial_n=[];
fixation_break_n=[];
artifact_PretargetToTarget_n=[];
artifact_BLintTo100msPostResponse_n=[];
artifact_BLintTo900ms_n=[];
artifact_neg1000_to_0ms_n=[];
hand_used = [];

erp = [];erp_HPF = []; erp_CSD=[]; Alpha = []; Pupil=[]; GAZE_X=[]; GAZE_Y=[]; n=0; artifacts_anywhereInEpoch = 0;
allRT=[]; allrespLR=[]; allTrig=[]; numtr=0; falsealarm=[];   % note allRT will be in sample points
    artifchansthistrial = [];
for f=1:length(files)
    disp(f)
    pathstemp = [path_temp filesep 'ProcessedSubjects/HAPPE' filesep dyscAMRM filesep 'processed' filesep files{f}]; % this data is high pass filtered using 1Hz
    % load original
    if HPFused==1000
        load(pathstemp)
        clear pathstemp
        EEG.data  = double(EEG.data );
    else
        EEG = pop_loadbv(paths{f},filesorig{f});
        loadbvSK_DN
        %%
        EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
        chan_locations = 'actiCAP65_ThetaPhi.elp';
        EEG=pop_chanedit(EEG, 'load',{chan_locations 'filetype' 'besa'});
        EEG = eeg_checkset( EEG );
        EEG = pop_eegfiltnew(EEG, [],0.8,[],0,[],0); %low pass filter of 0.8Hz, the third last zero says wherther high pass or low pass is used
        EEG.setname='lowpassEEG';
        EEG = eeg_checkset( EEG );
        EEGorig=EEG;
        % interpolate bad channels % in case any other bad channels are found
        %     % through thresholding
        %     if ~isempty(badchans)
        %         EEG.chanlocs = chanlocs;
        %         EEG=eeg_interp(EEG,[badchans],'spherical');
        %     end
        
        clear EEG
        load(pathstemp)
        clear pathstemp pathsorig
        EEG_HPF = EEG;
        y1temp = fft(EEG_HPF.data);
        y2temp = fft(EEGorig.data);
        
        EEG_HPF.data = double(ifft(y1temp+y2temp)); % add the low pass filtered data to the processed data
        EEG = EEG_HPF;
    end
    %% Ger's method to filter alpha:
    EEG_alpha.data = filtfilt(H1,G1,EEG.data')';
    
    %% Dan's method to filter alpha
    %%%%%%%%%%% band-pass filter, isolating the alpha band at 8-14 Hz. ftrans,
    %%%%%%%%%%% ftype etc are specialised filter settings for alpha
    %     EEG_alpha = pop_firpm(EEG, 'fcutoff', [bandlimits(1,1) bandlimits(1,2)], 'ftrans', 1, 'ftype', ...
    %         'bandpass', 'wtpass', 1, 'wtstop', 28.7822, 'forder', 1822); % Parks-McClellan filter
    %
    if HPF
        %EEG_HPF.data = double(EEG_HPF.data(1:nchan,:));
         %[H_HP1,G_HP1]=butter(4,2*HPFcutoff/old_fs,'high');   %Ger's old HPF method
         %EEG_HPF.data = filtfilt(H_HP1,G_HP1,EEG_HPF.data')'; %Ger's old HPF method
        %EEG_HPF = pop_eegfiltnew(EEG_HPF, 1,249,[],0,[],0);
        EEG.data = filtfilthd(HD_HPF,EEG.data')'; %New HPF method using filters designed with "fdatool" tool in MATLAB
        %EEG_HPF.data = eegfilt(EEG_HPF.data,old_fs,HPFcutoff,0);   
        disp('HPF finished')
    end
    
    % First LP Filter
  %  if LPF, EEG.data = eegfilt(EEG.data,old_fs,0,LPFcutoff); end 
    if LPF
    [H_LP1,G_LP1]=butter(4,2*LPFcutoff/old_fs,'low');   %Ger's old HPF method
    EEG.data = filtfilt(H_LP1,G_LP1,EEG.data')'; %Ger's old HPF method
    end

   % if LPF, EEG_HPF.data = eegfilt(EEG_HPF.data,old_fs,0,LPFcutoff); end    
    % average-reference the whole continuous data (safe to do this now after interpolation):
    %EEG.data = EEG.data - repmat(mean(EEG.data([1:nchan],:),1),[nchan,1]);
    %EEG_HPF.data = EEG_HPF.data - repmat(mean(EEG_HPF.data([1:nchan],:),1),[nchan,1]);
    %EEG_alpha.data = EEG_alpha.data - repmat(mean(EEG_alpha.data([1:nchan],:),1),[nchan,1]);
    
    
    %% Sync up the Eyelink data:
    if ~exist(ET_matfiles{f}, 'file') %DN: if ET matfile has NOT has been saved previouslty,
        FixEyelinkMessages %then calculate and save it now
    end
    load(ET_matfiles{f}) %DN: load the ET mat file, including eyelinkevent
    %Add an extra 4 rows into the EEG struct - 'TIME'
    %'GAZE_X' 'GAZE_Y' 'AREA'. This will add these as extra channels onto EEG.data
    %So the final channel is the pupil area (i.e. diameter):
    
    EEG = pop_importeyetracker(EEG,ET_matfiles{f},[first_event last_event]...
        ,[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},0,1,0,1);

    % detrend 
    %%
    numev = length(EEG.event);
    
    % Fish out the event triggers and times
    clear trigs stimes RT motion_on
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    % remove all trigs >106 since there should be only 6 triggers
    % anyway
    trigs = trigs(trigs<107);
    targtrigs = [];starttrig=[];
    for i=1:length(trigs)
        if any( trigs(i)==targcodes)
            targtrigs = [targtrigs i];
        end
        if any( trigs(i)==5)
            starttrig = [starttrig i];
        end        
    end
    %         if trigs(targtrigs(end))==trialCond(1)
    %             motion_on = targtrigs(1:end-1); % GL: indices of trigs when motion on. get rid of last trig, it was a repeat
    %         else
    motion_on = targtrigs;
    %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import eye tracking
    ETs = pop_importeyetracker(EEG,ET_matfiles{f},[first_event last_event]...
        ,[1:4] ,{'TIME' 'GAZE_X' 'GAZE_Y' 'AREA'},1,1,0,1);    
    Pupil_ch=length(EEG.data(:,1)); %Now that the eyelink data is added to the EEG struct Find the channel number for pupil area/diameter
    GAZE_X_ch=length(EEG.data(:,1))-2;
    GAZE_Y_ch=length(EEG.data(:,1))-1;
    Pupil_t = length(EEG.data(:,1))-3;
    
    ET_data     = ETs.data([Pupil_t; GAZE_X_ch; GAZE_Y_ch; Pupil_ch],:);
    ETevent  = ETs.event;
    
    %[ET_data] = interpolateET(ET_data, [stimes(motion_on)' trigs(motion_on)'] , ETevent, eyelinkevent, fs, ET_matfiles{f}, dataset);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:length(motion_on)-1 %don't count the last trial
        numtr = numtr+1;
        if n==1
            repeat01(numtr)=2;
        end        
        clear ep ep_CSD ep_HPF ep_alpha ep_art_reject ep_test ep_filt_Alpha_Hz ep_filt_abs_cut
        if ~isnan(motion_on(n))
            
            locktime = stimes(motion_on(n)); % Lock the epoch to coherent motion onset.
            try
                if motion_on(n)<length(trigs)
                    % false alarm: whenever they press between the motions
                    if ~isempty(trigs(starttrig(n)+1:motion_on(n)-1))
                        temp = trigs(starttrig(n)+1:motion_on(n)-1);
                        falsealarm(numtr) = length(temp(temp==12))+length(temp(temp==13));
                        clear temp;
                    else
                        falsealarm(numtr)=0;
                    end
                    % classify as repeat or no repeat
                    if n>1
                        m(1) = ismember(trigs(motion_on(n)),targcodesLR(:,1));
                        m(2) =  ismember(trigs(motion_on(n-1)),targcodesLR(:,1));
                        if  (m(2)-m(1))==0% is it the same side?
                            repeat01(numtr)=1; %not same side since =0
                        else
                            repeat01(numtr)=0;
                        end
                    end                
                    if trigs(motion_on(n)+1)==12 ||trigs(motion_on(n)+1)==13 %either left or right
                        response_time = stimes(motion_on(n)+1)-locktime; % time in samples from beginning of motion to response.
                        response_time = floor(response_time);
                        if response_time>default_response_time*fs
                            response_time = default_response_time*fs;
                        end
                    else
                        response_time = default_response_time*fs;
                    end
                else
                    response_time = default_response_time*fs;
                end
            end
            try
                ep = EEG.data(1:nchan,locktime+ts);   % chop out an epoch
                ep_CSD = CSD(ep, G_CSD, H_CSD); % CSD epoch
                ep_HPF = EEG.data(1:nchan,locktime+ts);
                ep_alpha = EEG_alpha.data(1:nchan,locktime+ts);

            catch
                disp('EEG ended too soon')
                allTrig(numtr) = 0;
                allrespLR(numtr) = 0;
                allRT(numtr) = 0;
                erp(:,:,numtr) = zeros(nchan,ERP_samps);
                erp_HPF(:,:,numtr) = zeros(nchan,ERP_samps);
                Alpha(:,:,numtr) = zeros(nchan,Alpha_samps);
                rejected_trial_n(numtr)=1;
                if ~strcmp(subject_folder(s),'TTN49')%for TTN49 we started EEG recording slightly too late so missed the first couple of EEG samples
                    if n ~= length(motion_on) %ignore it if this happens on final motion trigger, it was an extra/repeat and cut off
                        keyboard
                    end
                end
                continue;
            end
            try
                %     (1,:,:) = time
                %     (2,:,:) = X motion
                %     (3,:,:) = Y motion
                %     (4,:,:) = raw pupil diameter, AREA
                %     (5,:,:) = blink interpolated pupil diameter
                %     (6,:,:) = interpolated and low-pass filtered pupil diameter < 4 Hz
                %     (7,:,:) = interpolated and band-pass filtered pupil diameter 0.05 - 4
                %     (8,:,:) = instantaneous phase
                %     (9,:,:) = interpolated and low-pass filtered pupil diameter < 1 Hz
                %     (10,:,:) = interpolated and band-pass filtered pupil diameter 0.05 - 1 Hz
                % channels are [Pupil_t; GAZE_X_ch; GAZE_Y_ch; Pupil_ch]
                ep_pupil = ET_data(:,locktime+ts_pupil);   % chop out an epoch of pupil diameter
                ep_GAZE_X = ET_data(2,locktime+ts_pupil);
                ep_GAZE_Y = ET_data(3,locktime+ts_pupil);
                
            catch
                disp('Pupil Diameter data ended too soon')
                allTrig(numtr) = 0;
                allrespLR(numtr) = 0;
                allRT(numtr) = 0;
                Pupil(:,numtr) = zeros(ERP_samps);
                GAZE_X(:,numtr) = zeros(ERP_samps);
                GAZE_Y(:,numtr) = zeros(ERP_samps);
                rejected_trial_n(numtr)=1;
                keyboard
                continue;
            end
            
            % detrend
            BLamp = mean(ep_pupil(:,find(t_pupil>BLint(1) & t_pupil<BLint(2))),2);
            ep_pupil = ep_pupil - repmat(BLamp,[1,length(t_pupil)]); % baseline correction

            BLamp = mean(ep_GAZE_X(:,find(t_pupil>BLint(1) & t_pupil<BLint(2))),2);
            ep_GAZE_X = ep_GAZE_X - repmat(BLamp,[1,length(t_pupil)]); % baseline correction
            
            BLamp = mean(ep_GAZE_Y(:,find(t_pupil>BLint(1) & t_pupil<BLint(2))),2);
            ep_GAZE_Y = ep_GAZE_Y - repmat(BLamp,[1,length(t_pupil)]); % baseline correction            
            
            BLamp =mean(ep(:,find(t>BLint(1) & t<BLint(2))),2); % record baseline amplitude (t<0) for each channel,
            ep = ep - repmat(BLamp,[1,length(t)]); % baseline correction
            
            BLamp =mean(ep_CSD(:,find(t>BLint(1) & t<BLint(2))),2); % record baseline amplitude (t<0) for each channel,
            ep_CSD = ep_CSD - repmat(BLamp,[1,length(t)]); % baseline correction
            
            BLamp = mean(ep_HPF(:,find(t>BLint(1) & t<BLint(2))),2);
            ep_HPF = ep_HPF - repmat(BLamp,[1,length(t)]); % baseline correction
            
            BLamp = mean(ep_alpha(:,find(t>BLint(1) & t<BLint(2))),2);
            ep_alpha = ep_alpha - repmat(BLamp,[1,length(t)]); % baseline correction
            
            
            clear BLamp
            ep_test = [find(ts==-0.8*fs):find(ts==(0*fs))];
            if isempty(ep_test)
                disp('Empty epoch for art rejection')
                keyboard
            end
            ep_test = [find(t>BLint(1) & t<floor(((response_time*1000/fs)+100)))];
            if isempty(ep_test)
                disp('Empty epoch for art rejection2')
                keyboard
            end
            
            artifchans_thistrial = ARchans(find(max(abs(ep_HPF(ARchans,find(t<0))),[],2)>artifth | max(abs(ep_HPF(ARchans,find(t>BLint(1) & t<floor(((response_time*1000/fs)+100))))),[],2)>artifth));
            
            artifchans_blinks_thistrial = ARchans(find(max(abs(ep_HPF(ARchans_for_blinks,find(t<0))),[],2)>artifth));
            
            artifchans_blinks_thistrial(find(ismember(artifchans_blinks_thistrial,artifchans_thistrial))) = [];

            artifchans_thistrial = [artifchans_thistrial,artifchans_blinks_thistrial(find(~ismember(artifchans_blinks_thistrial,artifchans_thistrial)))];
            
            artifchans_PretargetToTarget_thistrial = ARchans(find(max(abs(ep_HPF(ARchans,find(ts==PretargetARwindow(1)*fs):find(ts==0))),[],2)>artifth));  % pre-target artifact rejection from -500-0ms only [find(ts==-.500*fs) gives you the point in samples -500ms before the target]
            artifchans_BLintTo100msPostResponse_thistrial = ARchans(find(max(abs(ep_HPF(ARchans,find(t>BLint(1) & t<floor(((response_time*1000/fs)+500))))),[],2)>artifth)); %Baseling until 100ms after response.
            
            %blinkChans = [1 2 3 7 33 36 41 46];
            %artifchans_BLintTo100msPostResponse_thistrial =...
            %    artifchans_BLintTo100msPostResponse_thistrial(~ismember(artifchans_BLintTo100msPostResponse_thistrial,blinkChans));
     
            artifchans = [artifchans artifchans_BLintTo100msPostResponse_thistrial];
            
            temp = zeros([1 65]); temp(artifchans_BLintTo100msPostResponse_thistrial)=1;
            artifchansthistrial = [artifchansthistrial; temp]; % for trials  
            clear temp
            
            artifchans_BLintTo900ms_thistrial = ARchans(find(max(abs(ep_HPF(ARchans,find(ts==-0.1*fs):find(ts==0.9*fs))),[],2)>artifth));  % artifact rejection from -100 to 900ms only
            
            artifchans_neg1000_to_0ms_thistrial= ARchans(find(max(abs(ep_HPF(ARchans,find(ts==-1*fs):find(ts==0))),[],2)>artifth));  % artifact rejection from -1000 to 0ms only
            

            % store the artifact in blocks and trials
            if ~isempty(artifchans_thistrial)
                artifacts_anywhereInEpoch = artifacts_anywhereInEpoch+1;
            end   % artifact rejection (threshold test)
            
            if artifchans_PretargetToTarget_thistrial
                artifact_PretargetToTarget_n(numtr)=1;
            else
                artifact_PretargetToTarget_n(numtr)=0;
            end
            
            if artifchans_BLintTo100msPostResponse_thistrial
                artifact_BLintTo100msPostResponse_n(numtr)=1;
            else
                artifact_BLintTo100msPostResponse_n(numtr)=0;
            end
            
            if artifchans_BLintTo900ms_thistrial
                artifact_BLintTo900ms_n(numtr)=1;
            else
                artifact_BLintTo900ms_n(numtr)=0;
            end
            
            if artifchans_neg1000_to_0ms_thistrial
                artifact_neg1000_to_0ms_n(numtr)=1;
            else
                artifact_neg1000_to_0ms_n(numtr)=0;
            end
            
            
            % scres = 1024 x 768: 512, 384 is middle. 3 deg is 76 pixels. Nope!
            % one pixel is 53/1280=0.04cm,  41cm wide and 30cm high, 5
            % degrees is 8cm = 200 pixels
            rangerx = 100; rangery = 76;
            artif_ET_pretarg = find(ep_pupil(2,find(ts<=0))<-rangerx | ep_pupil(2,find(ts<=0))>rangerx ...
                                   |ep_pupil(3,find(ts<=0))<-rangery | ep_pupil(3,find(ts<=0))>rangery);
                               
           pupilxtemp = ep_pupil(2,find(ts_pupil>=0*fs & ts_pupil<=response_time));
           pupilytemp = ep_pupil(3,find(ts_pupil>=0*fs & ts_pupil<=response_time));

           pupilxtemp = pupilxtemp-mean(pupilxtemp);
           pupilytemp = pupilytemp-mean(pupilytemp);
                              
%            
            artif_ET_BL_resp = find(abs(pupilytemp)>rangery);
            %pupilAll(:,:,numtr) = [pupilxtemp;pupilytemp];  clear pupilxtemp pupilytemp
            % 0 = reject, 1 = keep
            if length(artif_ET_pretarg) > 0, ET_pretarg_artrej(numtr) = 0; else ET_pretarg_artrej(numtr) = 1; end
            if length(artif_ET_BL_resp) > 0, ET_BL_resp_artrej(numtr) = 0; else ET_BL_resp_artrej(numtr) = 1; end
            
            ep = double(ep); % filtfilt needs doubles
            ep_HPF = double(ep_HPF);
            %%  Ger's method for alpha::
            ep_filt_Alpha_Hz = filtfilt(H1,G1,ep')'; % alpha filter again.
            %% Dan's Method:
            %             ep_filt_Alpha_Hz =double(ep_alpha);
            %%
            ep_pupil = double(ep_pupil);
            ep_GAZE_X = double(ep_GAZE_X);
            ep_GAZE_Y = double(ep_GAZE_Y);
            
            %%%%%%%% rectifying the data and chopping off ends %%%%%%%%
            ep_filt_abs_cut = abs(ep_filt_Alpha_Hz(:,find(ts==ts_crop(1)):find(ts==ts_crop(end)))); % 64x701
            % Smoothing. This goes from 1:700, leaving out the final
            % sample, 0ms.
            alpha_temp = []; Alpha_smooth_time = []; Alpha_smooth_sample = [];
            for q = 1:size(ep_filt_abs_cut,1)
                counter = 1;
                for windowlock = 26:25:size(ep_filt_abs_cut,2)-25 % 1,26,51,etc. 26 boundaries = 1:50, 51 boundaries = 26:75, 676 boundaries = 651:
                    alpha_temp(q,counter) = mean(ep_filt_abs_cut(q,windowlock-25:windowlock+24));
                    Alpha_smooth_time(counter) = t_crop(windowlock);
                    Alpha_smooth_sample(counter) = ts_crop(windowlock);
                    counter = counter+1;
                end
            end
            
            %             figure
            %             plot(ep_filt_abs_cut(54,:))
            %             figure
            %             plot(Alpha_smooth_time,alpha_temp(54,:))
            %             keyboard
            
            %             figure, hold on, for i = 1:64, plot(ep(i,find(ts==ts_crop(1)):find(ts==ts_crop(end))),'b'), end; keyboard
            erp(:,:,numtr) = ep(:,find(ts==ts_crop(1)):find(ts==ts_crop(end)));
            erp_CSD(:,:,numtr) = ep_CSD(:,find(ts==ts_crop(1)):find(ts==ts_crop(end))); %For CSD
            erp_HPF(:,:,numtr) = ep_HPF(:,find(ts==ts_crop(1)):find(ts==ts_crop(end)));
            
            Pupil(:,:,numtr)= ep_pupil(:, find(ts_pupil==ts_pupilcrop(1)):find(ts_pupil==ts_pupilcrop(end)));
            GAZE_X(:,numtr)= ep_GAZE_X(find(ts_pupil==ts_pupilcrop(1)):find(ts_pupil==ts_pupilcrop(end)));
            GAZE_Y(:,numtr)= ep_GAZE_Y(find(ts_pupil==ts_pupilcrop(1)):find(ts_pupil==ts_pupilcrop(end)));
            
            Alpha(:,:,numtr) = alpha_temp;
            allTrig(numtr) = trigs(motion_on(n));
            
            load(matfiles{f},'par');
            
            try % change this
                rejected_trial_n(numtr)=0;
                if trigs(motion_on(n)+1)==12 && par.responsehand==1 % they pressed the right button
                    allrespLR(numtr) = 1;
                    allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
                    hand_used(numtr) = 1;
                elseif trigs(motion_on(n)+1)==13 && par.responsehand==2% they pressed the right button
                    allrespLR(numtr) = 1;
                    allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));
                    hand_used(numtr) = 2;
                elseif trigs(motion_on(n)+1)==12 && par.responsehand==2% they pressed the wrong button
                    allrespLR(numtr) = 2;
                    allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));  
                    hand_used(numtr) = 1;
                elseif trigs(motion_on(n)+1)==13 && par.responsehand==1% they pressed the wrong button
                    allrespLR(numtr) = 2;
                    allRT(numtr) = stimes(motion_on(n)+1)-stimes(motion_on(n));  
                    hand_used(numtr) = 2;
                else
                    allrespLR(numtr) = 3; % no response, to mark it out from artifact trials.
                    allRT(numtr) = 0;
                    hand_used(numtr)=0; %no hand used
                end
            catch
                rejected_trial_n(numtr)=1;
                allrespLR(numtr) = 0;
                allRT(numtr) = 0;
            end
            
        else
            rejected_trial_n(numtr)=1; %if isnan(motion_on(numtr)) - this should never happen since motion_on will always be a number
            allTrig(numtr) = 0;
            allrespLR(numtr) = 0;
            allRT(numtr) = 0;
            erp(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_HPF(:,:,numtr) = zeros(nchan,ERP_samps);
            Alpha(:,:,numtr) = zeros(nchan,Alpha_samps);
            Pupil(:,numtr) = zeros(1,ERP_samps);
            GAZE_X(:,numtr) = zeros(1,ERP_samps);
            GAZE_Y(:,numtr) = zeros(1,ERP_samps);
        end
    end
end


Alpha = Alpha(1:nchan,:,:);
erp = erp(1:nchan,:,:);
erp_HPF = erp_HPF(1:nchan,:,:);
Pupil=Pupil(:,:,:);
GAZE_X=GAZE_X(:,:);
GAZE_Y=GAZE_Y(:,:);
figure;
hist(artifchans,[1:nchan]); title([allsubj{s} ': ' num2str(artifacts_anywhereInEpoch) ' artifacts = ',num2str(round(100*(artifacts_anywhereInEpoch/length(allRT)))),'%']) % s from runafew
disp([allsubj{s},' number of trials: ',num2str(length(allRT))])
if length(allRT)~=size(Alpha,3)
    disp(['WTF ',allsubj{s},' number of trials: ',num2str(length(allRT)),' not same as Alpha'])
    keyboard
end


baseline_erp = mean(erp(:,find(t>=BLint(1) & t<=BLint(2)),:),2);
erp_temp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp

erp_temp = squeeze(mean(erp_temp(:,:,find(artifact_BLintTo100msPostResponse_n==0)),3));
%h1 = figure
%plottopo(erp_temp(:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
%    min(min(erp_temp(:,:)))  100], ...
%    'title',['ERP'],'ydir',1)
%saveas(h1, ['Figures/capplot_', num2str(subj) '.png'])
h1 = figure
plot(sum(artifchansthistrial'))
xticks(1:18:216)
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12'})
 title([subject_folder{1}, 'bad channels over trials' num2str(round(sum(ET_BL_resp_artrej & ~artifact_BLintTo100msPostResponse_n)/length(ET_BL_resp_artrej)*100))])
saveas(h1, ['Figures/badchans_', num2str(subj) '.png'])
[counts,centers] = hist(artifchans,[1:nchan]);
h1 = figure;
topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers','maplimits', [0 10]);
title([subject_folder{1},'no of trials, blinking trials ' num2str(sum(~ET_BL_resp_artrej))])
pause(1)
saveas(h1, ['Figures/caprej', num2str(subj) '.png'])

if HPFused==1000
    pathsave = [path_temp 'ProcessedSubjects/HAPPE/' dyscAMRM '/' allsubj{s}, '_' num2str(LPFused) 'fftalldots']
else
    pathsave = [path_temp 'ProcessedSubjects/' num2str(HPFused) 'Hz/' dyscAMRM '/' allsubj{s}, '_' num2str(LPFused) 'fftalldots'];
end
save(pathsave,...
    'Alpha','erp','erp_CSD','erp_HPF','Pupil','GAZE_X','GAZE_Y','allRT','allrespLR','allTrig','falsealarm', ...
    'artifchans','t_crop','Alpha_smooth_time','Alpha_smooth_sample',...
    'artifact_PretargetToTarget_n','artifact_BLintTo100msPostResponse_n',...
    'fixation_break_n','rejected_trial_n','artifact_BLintTo900ms_n','artifact_neg1000_to_0ms_n','hand_used','repeat01',...
    'ET_BL_resp_artrej', 'ET_pretarg_artrej')
close all
return;


% %% STFT:
% nchan=65;
% fs=500;
% clear stftC
% TC = [-1500:100:700];%100ms sliding window
% fftlen = 300;
% F = [0:20]*fs/fftlen; %Frequencies
% for tt=1:length(TC)
%     temp = abs(fft(erp(:,find(t_crop>=TC(tt),1)-fftlen/2+[1:fftlen],:),[],2))./(fftlen/2);
%     stftC(:,:,tt,:) = reshape(temp(:,1:length(F),:),[nchan length(F) 1 size(erp,3)]);
% end %stftC(Electrode,Frequency,Time,Trial)
%
% %Isolate time-range and collapse accross it
% Trange = find(TC>0 & TC<800);
% spec =squeeze(mean(stftC(:,:,Trange,:),3)); %spec(Electrode,Frequency,Trial)
% %Isolate frequency band within that time range and collapse accross it:
% band = find(F>8 & F<14);
%
% PreAlpha=squeeze(mean(spec(:,band,:),2)); %PreAlpha(electrode trial)
%
%
% Alpha_simon=squeeze(mean(stftC(:,band,:,:),2)); %Alpha_simon(Electrode,Time,Trial)
% %






