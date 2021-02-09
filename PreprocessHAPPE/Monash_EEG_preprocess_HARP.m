%% this is a custom script using the HAPPE preprocessing method which is designed for children and any person who are unable to fixate
fs=500;
targcodes = [101:106];

old_fs = 500; % old sample ratenameChanges
targcodesLR = zeros(3,2);
targcodesLR(:,1) = [101,103,105];
targcodesLR(:,2) = [102,104,106];

% in sample points, the ERP epoch
ts = -1.2*fs:2.2*fs; % -1000ms to 1880ms with 200ms either side.
t = ts*1000/fs;
ts_crop = -1*fs:2*fs;
t_crop = ts_crop*1000/fs;
BLint = [-100 0];   % baseline interval in ms
default_response_time = 1.1-0.1;

Alpha_samps = ((abs(t_crop(1))+abs(t_crop(end)))/50)-1; %50ms alpha samples
ERP_samps = length(t_crop);

nchan = 65;

% automatically high passes with 1Hz
LPFcutoff=20;       % Low Pass Filter cutoff
HPFcutoff=0.01;       % High Pass Filter - either 0.01, 0.1, or 0.25 for cuttoff as I have filters .mats designed for these using "fdatool" tool in MATLAB
LPF = 1;    % 1 = low-pass filter the data, 0=don't.

bandlimits(1,1) = 8; % defining the filter for alpha bandpass.
bandlimits(1,2) = 13;

[H1,G1]=butter(4,[2*(bandlimits(1,1)/old_fs) 2*(bandlimits(1,2)/old_fs)]); % alpha bandpass
[H2,G2]=butter(4,[2*(bandlimits(1,1)/fs) 2*(bandlimits(1,2)/fs)]); % alpha bandpass

PretargetARwindow=[-0.500,0];%time window (in seconds, must be factor of the 50ms alpha samples) to search for pre-target artifacts

% 1. enter path to the folder that has the datasets you want to analyze

src_folder_name=[path_temp 'ProcessedSubjects/' dyscAMRM '/'];


% 3. list channels of interest, including the 10-20 channels. User defined channels occur at the end of the sequence e.g. 'E39' 'E40'
%the 18 "10-20" channels that NEED to be in the chan_IDs: 'FP1' 'FP2' 'F3'
% 'F4' 'F7' 'F8' 'C3' 'C4' 'T3' 'T4' 'PZ' 'O1' 'O2' 'T5' 'T6' 'P3' 'P4' 'Fz'
chan_IDs={'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';...
    'T3';'C3';'Cz';'C4';'T4';'TP9';'CP5';'CP1';'CP2';'CP6';'TP10';'T5';...
    'P3';'Pz';'P4';'T6';'PO9';'O1';'Oz';'O2';'PO10';'AF7';'AF3';'AF4';...
    'AF8';'F5';'F1';'F2';'F6';'FT9';'FT7';'FC3';'FC4';'FT8';'FT10';'C5';'C1';'C2';'C6';...
    'TP7';'CP3';'CPz';'CP4';'TP8';'P5';'P1';'P2';'P6';'PO7';'PO3';'POz';'PO4';'PO8';'FCz'}';


% for MARA ICA decomposition
nameChanges.numbers = [12 16, 23, 27];
nameChanges.labels = {'T3','T4','T5','T6'};
% 4. run HAPPE in semi-automated setting with visualizations (=1) or fully-automated, no visualizations setting ( = 0)
pipeline_visualizations_semiautomated = 0;
%if semi-automated, what is the minimum and maximum frequency you want to visualize in the power spectrum figure for each file?
vis_freq_min = 2;
vis_freq_max = 57;
%if semi-automated, which frequencies do you want to generate spatial topoplots for within the figure?
freq_to_plot = [6 10 22 42 55];

% 5. for resting-state EEG, set task_EEG_processing = 0
%for task-related EEG, set task_EEG_processing = 1
task_EEG_processing = 1;

%if task-related EEG:
task_conditions = {'near', 'devr'}; %enter the stimulus condition tags

%if resting-state EEG:
% list all potential names of the matlab variable that contains the EEG data for your files:
% Note that variable names that include the string 'Segment' may cause
% difficulties with pop_importegimat (as data should not be pre-epoched).
% Consider renaming these variables if difficulties arise
potential_eeg_var_names = {'Category_1_Segment1','Category_1'};

% 6. do you want to segment your data? yes (=1) or no (=0)
segment_data = 0;

%if you are segmenting your task-related EEG:
%parameters to segment the data for each stimulus, in seconds:
task_segment_start = -0.5;
task_segment_end = 1.5;

%if you are segmenting your resting-state EEG:
%how long do you want your segments to be in seconds? here, 2 seconds is the default
segment_length = 2;

% 7. do you want to interpolate the specific channels' data determined to be artifact/bad within each segment?
%yes = 1, no = 0.
%This is segment-level channel interpolation from the FASTER EEGlab plug-in.
segment_interpolation = 0;

% 8. do you want to do segment rejection (using amplitude and joint probability criteria)?
%yes = 1, no = 0.
segment_rejection = 0;

% if you are rejecting segments, what minimum/maximum signal amplitude do you want to use as the artifact threshold?
reject_min_amp = -100;
reject_max_amp = 100;

% do you want to do segment rejection using all user-specified channels above ( = 0) or a subset of channels in an ROI ( = 1)?
ROI_channels_only = 0;

% if you want to do ROI segment rejection, which channels should be used in the ROI?
% ex ROI_channels = {'E27','E20'};
ROI_channels = {'E27','E20'};

% 9. Select the type of re-referencing you want. Average re-reference (=1)
% or re-referencing to another channel/subset of channels (=0)
average_rereference = 1;

% if you are referencing to another channel/subset of channels, what are they?
% make sure to use the channel name given above
% ex ROI_channels = {'E57','E100'};
NO_AVERAGE_REREF_channel_subset = {'E57','E100'};

% 10. Select the format to save your processed data at the end of HAPPE!
%save_as_format = 1 will save the processed data as a .txt file.(electrodes as columns, time as rows)
%save_as_format = 2 will save the processed data as a .mat file (matlab format)
%save_as_format = 3 will save the processed data as a .set file (EEGlab format)
save_as_format = 2;

%~~~~~~~~~~~~~~~~~~~~~~ no need to edit beyond this point ~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% make output folders:

%make folder for intermediate waveleted files
savepath{1}=[src_folder_name filesep 'intermediate1_wavclean' filesep  allsubj{s}];
savepath{2}=[src_folder_name filesep 'intermediate2_ICAclean' filesep  allsubj{s}];
savepath{3}=[src_folder_name filesep 'intermediate3_segmented' filesep  allsubj{s}];
savepath{4}=[src_folder_name filesep 'processed' filesep  allsubj{s}];
savepath{5}=[src_folder_name filesep 'processed' filesep  allsubj{s} filesep 'figures'];
if ~isdir (savepath{1})
    mkdir (savepath{1});
end

%make folder for post-ICA, uninterpolated files
if ~isdir (savepath{2})
    mkdir (savepath{2});
end

%make folder for segment-level files (if user selects segmentation option)
if ~isdir (savepath{3})
    mkdir (savepath{3});
end
%make folder for final preprocessed files
if ~isdir (savepath{4})
    mkdir (savepath{4});
end
%make folder for figures
if ~isdir (savepath{5})
    mkdir (savepath{5});
end

%% add relevant folders to path

% add HAPPE script path
happe_directory_path = 'happe-master';

% will eventually allow users to set own eeglab path -- for now, assume
% using eeglab14_0_0b included in HAPPE
eeglab_path = [happe_directory_path filesep 'Packages' filesep 'eeglab14_0_0b'];

% add HAPPE subfolders and EEGLAB plugin folders to path
addpath([happe_directory_path filesep 'acquisition_layout_information'],[happe_directory_path filesep 'scripts'],...
    eeglab_path,genpath([eeglab_path filesep 'functions']));
rmpath(genpath([eeglab_path filesep 'functions' filesep 'octavefunc']));

plugin_directories = dir([eeglab_path filesep 'plugins']);
plugin_directories = strcat(eeglab_path,filesep,'plugins',filesep,{plugin_directories.name},';');
addpath([plugin_directories{:}]);

% % add cleanline path
if exist('cleanline','file')
    cleanline_path = which('eegplugin_cleanline.m');
    cleanline_path = cleanline_path(1:findstr(cleanline_path,'eegplugin_cleanline.m')-1);
    addpath(genpath(cleanline_path));
else
    error('Please make sure cleanline is on your path');
end


% frontal channels, occipital channels, and POz and Pz, CPz, CP1, CP2, P1, P2

ARchans = [1:65];
ARchans_for_blinks = [1:65];
artifth = 100;
artifchans=[];  % keep track of channels on which the threshold is exceeded, causing trial rejection

% chanlocs = readlocs('cap128.loc');
chan_locations = 'actiCAP65_ThetaPhi.elp';
chanlocs=readlocs ('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP with reference channel 'FCz' included - hence 65 chans
%% go to the folder with data, pull file names to feed script,and initialize the arrays to store file specific data quality metrics
% cd get file list
FileNames = files;
% intialize report metrics
chan_index=[1:length(chan_IDs)];
Number_ICs_Rejected=[];
Number_Good_Channels_Selected=[];
Interpolated_Channel_IDs=[];
Percent_ICs_Rejected=[];
Percent_Variance_Kept_of_Post_Waveleted_Data=[];
File_Length_In_Secs=[];
Number_Channels_User_Selected=[];
Percent_Good_Channels_Selected=[];
Median_Artifact_Probability_of_Kept_ICs=[];
Mean_Artifact_Probability_of_Kept_ICs=[];
Range_Artifact_Probability_of_Kept_ICs=[];
Min_Artifact_Probability_of_Kept_ICs=[];
Max_Artifact_Probability_of_Kept_ICs=[];
Number_Segments_Post_Segment_Rejection=[];
src_file_ext = '.vhdr';
%iterate the following preprocessing pipeline over all your data files:
for current_file = 1:length(FileNames)
    processedfile = [savepath{4} filesep strrep(FileNames{current_file}, src_file_ext,'_processed.mat')];
    if exist(processedfile,'file')==2
        clear processedfile
        current_file
    else
        clear processedfile
        current_file
        %% load file and get sampling rate, save with double precision
        
        % import data into eeglab and store the file's length in seconds for outputting later
        if task_EEG_processing == 1
            EEG = pop_loadbv(paths{current_file},files{current_file});
            %% From loadbvSK that simon wrote
            loadbvSK_DN
            filepath = [paths{current_file} filesep FileNames{current_file}];
            %%
            EEGloaded = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
            clear EEG
            %EEGloaded = pop_readegi(filepath, [],[],'auto');
            events=EEGloaded.event;
            complete_event_info=EEGloaded.urevent;
            srate=double(EEGloaded.srate);
            
            
            % interpolate the channels which are 100% bad (i.e. turned off)
            if ~isempty(badchans)
                EEGloaded.chanlocs = chanlocs;
                EEGloaded=eeg_interp(EEGloaded,[badchans],'spherical');
            end
        elseif task_EEG_processing == 0
            load(FileNames{current_file});
            srate=double(samplingRate);
            file_eeg_vname = intersect(who,potential_eeg_var_names);
            
            try
                EEGloaded = pop_importegimat(FileNames{current_file}, srate, 0.00, file_eeg_vname{1});
            catch err_msg
                if strcmp(err_msg.identifier,'MATLAB:badsubscript')
                    error('sorry, could not read the variable name of the EEG data, please check your file')
                else
                    error(err_msg.message);
                end
            end
        end
        
        EEGloaded.setname='rawEEG';
        EEG = eeg_checkset(EEGloaded);
        File_Length_In_Secs(current_file)=EEG.xmax;
        
        % edit channel locations (does not import properly from netstation by default)
        EEG=pop_chanedit(EEG, 'load',{chan_locations 'filetype' 'besa'});
        EEG = eeg_checkset( EEG );
        
        %load 10-20 EEG system labels for electrode names (for MARA to reference)
        for i=1:length(nameChanges.numbers)
            EEG=pop_chanedit(EEG, 'changefield',{nameChanges.numbers(i)  'labels' nameChanges.labels{i}});
        end
        EEG = eeg_checkset( EEG );
        
        %% filter the data with 1hz highpass (for srate 250), bandpass 1hz-249hz (for srate 500, ICA doesn't reliably work well with frequencies above 250hz)
        if srate<500
            EEG = pop_eegfiltnew(EEG, [],1,[],1,[],0);
        elseif srate >= 500
            EEG = pop_eegfiltnew(EEG, 1,249,[],0,[],0);
        end
        EEG.setname='rawEEG_f';
        EEG = eeg_checkset( EEG );
        
        %% select EEG channels of interest for analyses and 10-20 channels from the list you specified at the top of the script
        
        EEG = pop_select( EEG,'channel', chan_IDs);
        EEG.setname='rawEEG_f_cs';
        EEG = eeg_checkset( EEG );
        full_selected_channels = EEG.chanlocs;
        
        
        %     %% Define channels, having combined Brain Products and Biosemi data
        %     chanlocs = EEG.chanlocs;
        %     plot_chans = 1:65;
        %     % exclude_chans = [];
        %     %plot_chans = [1:64];
        %     left_hemi = [1,33,34,3,37,4,38,9,43,8,42,41,12,47,13,48,19,52,18,51,17,23,56,24,57,61,60,28,29];
        %     right_hemi = [2,35,36,39,6,40,7,10,44,11,45,46,49,15,50,16,20,54,21,55,22,58,26,59,27,63,64,51,32];
        %     centre_chans = [5,65,14,53,25,62,30];
        %     elec_pairs = [1,2;33,36;34,35;3,7;37,40;4,6;38,39;...
        %         41,46;42,45;8,11;43,44;9,10;48,49;13,15;47,50;12,16;
        %         17,22;51,55;18,21;52,54;19,20;...
        %         23,27;56,59;24,26;57,58;60,64;61,63;29,31;28,32];
        %
        %     tester = zeros(65,1);
        %     figure
        %     topoplot(tester,chanlocs,'maplimits', ...
        %         [min(tester)  max(tester)],'electrodes','numbers','plotchans',plot_chans);
        %     figure
        %     topoplot(tester,chanlocs,'maplimits', ...
        %         [min(tester)  max(tester)],'electrodes','labels','plotchans',plot_chans);
        %% reduce line noise in the data (note: may not completely eliminate, re-referencing helps at the end as well)
        EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',chan_index,'computepower',1,'linefreqs',...
            [60 120] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
            'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
        EEG.setname='rawEEG_f_cs_ln';
        EEG = eeg_checkset(EEG);
        
        % close window if visualizations are turned off
        if pipeline_visualizations_semiautomated == 0
            close all;
        end
        %% crude bad channel detection using spectrum criteria and 3SDeviations as channel outlier threshold, done twice
        
        EEG = pop_select(EEG,'nochannel',65);
        EEG = pop_rejchan(EEG, 'elec',chan_index(1:64),'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]); %65 is zero and will be ommited from the ICA
        EEG.setname='rawEEG_f_cs_ln_badc';
        EEG = eeg_checkset( EEG );
        
        EEG = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],'norm','on','measure','spec','freqrange',[1 125]);
        EEG.setname='rawEEG_f_cs_ln_badc2x';
        EEG = eeg_checkset( EEG );
        selected_channel_locations=EEG.chanlocs;
        
        %save the names of the rejected channels for output table after the pipeline finishes
        selected_channel_labels={selected_channel_locations.labels};
        bad_channels_removed= setdiff(chan_IDs, selected_channel_labels);
        [~,ROI_indices_in_selected_chanlocs] = intersect(selected_channel_labels,ROI_channels);
        
        %% run wavelet-ICA (ICA first for clustering the data, then wavelet thresholding on the ICs)
        %uses a soft, global threshold for the wavelets, wavelet family is coiflet (level 5), threshold multiplier .75 to remove more high frequency noise
        %for details, see wICA.m function
        
        try
            if pipeline_visualizations_semiautomated == 0
                [wIC, A, W, IC] = wICA(EEG,'runica', 1, 0, srate, 5);
            elseif pipeline_visualizations_semiautomated == 1
                [wIC, A, W, IC] = wICA(EEG,'runica', 1, 1, srate, 5);
            end
        catch wica_err
            if strcmp ('Output argument "wIC" (and maybe others) not assigned during call to "wICA".',wica_err.message)
                error('Error during wICA, most likely due to memory settings. Please confirm your EEGLAB memory settings are set according to the description in the HAPPE ReadMe')
            else
                rethrow(wica_err)
            end
        end
        
        %reconstruct artifact signal as channelsxsamples format from the wavelet coefficients
        artifacts = A*wIC;
        
        %reshape EEG signal from EEGlab format to channelsxsamples format
        EEG2D=reshape(EEG.data, size(EEG.data,1), []);
        
        %subtract out wavelet artifact signal from EEG signal
        wavcleanEEG=EEG2D-artifacts;
        
        %save wavelet cleaned EEG data file to folder with extension _wavclean.mat
        
        save([savepath{1} filesep strrep(FileNames{current_file}, src_file_ext,'_wavclean.mat')],'wavcleanEEG')
        save([savepath{1} filesep strrep(FileNames{current_file}, src_file_ext,'_prewav.mat')],'EEG2D')
        
        %% reimport into EEGlab
        EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[savepath{1} filesep strrep(FileNames{current_file},src_file_ext,'_wavclean.mat')],'srate',srate,'pnts',0,'xmin',0,'chanlocs',selected_channel_locations);
        EEG.setname='wavcleanedEEG';
        EEG = eeg_checkset( EEG );
        
        % import event tags if needed
        if task_EEG_processing == 1
            EEG.event=events;
            EEG.urevent=complete_event_info;
        end
        
        %% run ICA to evaluate components this time
        ranknumbs = rank(EEG.data);
        EEG = pop_runica(EEG, 'extended',1,'interupt','on');
        clear n;
        EEG = eeg_checkset( EEG );
        %save the ICA decomposition intermediate file before cleaning with MARA
        EEG = pop_saveset(EEG, 'filename',strrep(FileNames{current_file}, src_file_ext,'_ICA.set'),'filepath',[savepath{2}]);
        clear EEG
        %% use MARA to flag artifactual IComponents automatically if artifact probability > .5
        EEG = pop_loadset(strrep(FileNames{current_file}, src_file_ext,'_ICA.set'),savepath{2});
        [~,EEG,~]=processMARA ( EEG,EEG,EEG, [0, 0, pipeline_visualizations_semiautomated,...
            pipeline_visualizations_semiautomated , pipeline_visualizations_semiautomated] );
        
        EEG.reject.gcompreject = zeros(size(EEG.reject.gcompreject));
        EEG.reject.gcompreject(EEG.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
        EEG.setname='wavcleanedEEG_ICA_MARA';
        EEG = eeg_checkset( EEG );
        
        % store MARA related variables to assess ICA/data quality
        index_ICs_kept=(EEG.reject.MARAinfo.posterior_artefactprob < 0.5);
        median_artif_prob_good_ICs = median(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        mean_artif_prob_good_ICs = mean(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        range_artif_prob_good_ICs = range(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        min_artif_prob_good_ICs = min(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        max_artif_prob_good_ICs = max(EEG.reject.MARAinfo.posterior_artefactprob(index_ICs_kept));
        for i=1:EEG.trials
            EEG.icaact(:,:,i)=(EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:,i);
        end
        %store IC variables and calculate variance of data that will be kept after IC rejection:
        ICs_to_keep =find(EEG.reject.gcompreject == 0);
        ICA_act = EEG.icaact;
        ICA_winv =EEG.icawinv;
        
        %variance of wavelet-cleaned data to be kept = varianceWav:
        [projWav, varianceWav] =compvar(EEG.data, ICA_act, ICA_winv, ICs_to_keep);
        
        %% reject the ICs that MARA flagged as artifact
        artifact_ICs=find(EEG.reject.gcompreject == 1);
        EEG = pop_subcomp( EEG, artifact_ICs, 0);
        EEG.setname='wavcleanedEEG_ICA_MARA_rej';
        EEG = eeg_checkset( EEG );
        
        %% save the post-MARA cleaned intermediate file before interpolating anything
        EEG = pop_saveset(EEG, 'filename',strrep(FileNames{current_file}, src_file_ext,'_ICAcleanedwithMARA.set'),'filepath',savepath{2});
        clear EEG
        %% segment data according to data type
        % to add Ger's code here
        EEG = pop_loadset(strrep(FileNames{current_file}, src_file_ext,'_ICAcleanedwithMARA.set'),savepath{2});
        allEvent = cell2mat({EEG.event.type}');
        targcodes = {'101','102','103','104','105','106'};
        idxtrigger=find(allEvent>100);
        
        [EEG.event(idxtrigger(end)).type] = 10;
        
        
        if segment_data
            EEG = pop_epoch(EEG, targcodes,[-0.5 2.0], 'verbose', 'no', 'epochinfo', 'yes');
        end
        %     %%
        %     if segment_data
        %         if ~task_EEG_processing
        %             EEG=eeg_regepochs(EEG,'recurrence',segment_length,'limits',[0 segment_length], 'rmbase', [NaN]);
        %         else
        %             EEG = pop_epoch(EEG, task_conditions, [task_segment_start task_segment_end], 'verbose', 'no', 'epochinfo', 'yes');
        %         end
        %     end
        
        EEG = pop_saveset( EEG, 'filename',strrep(FileNames{current_file},src_file_ext,'_segmented.set'),'filepath',savepath{3});
        %% if selected option, interpolate bad data within segments from "good channels" only:
        if segment_interpolation
            
            %use only the good channels to evaluate data:
            eeg_chans=[1:length(selected_channel_locations)];
            
            %evaluate the channels for each segment and interpolate channels with bad
            %data for that each segment using the FASTER program, interpolating channels scoring above/below z threshold of 3 for an segment:
            ext_chans=[];
            o.epoch_interp_options.rejection_options.measure = [1 1 1 1];
            o.epoch_interp_options.rejection_options.z = [3 3 3 3];
            
            if  length(size(EEG.data)) > 2
                status = '';
                lengths_ep=cell(1,size(EEG.data,3));
                for v=1:size(EEG.data,3)
                    list_properties = single_epoch_channel_properties(EEG,v,eeg_chans);
                    lengths_ep{v}=eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options)));
                    status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) sprintf('\n')];
                end
                EEG=h_epoch_interp_spl(EEG,lengths_ep,ext_chans);
                EEG.saved='no';
                
                %add the info about which channels were interpolated for each segment to the EEG file
                EEG.etc.epoch_interp_info=[status];
            end
            
            EEG = pop_saveset(EEG, 'filename',strrep(FileNames{current_file}, src_file_ext,'_segments_interp.set'),'filepath',savepath{3});
        end
        
        %% rejection of bad segments using amplitude-based and joint probability artifact detection
        if segment_rejection
            if ROI_channels_only == 0
                EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan] ,[reject_min_amp],[reject_max_amp],[EEG.xmin],[EEG.xmax],2,0);
                EEG = pop_jointprob(EEG,1,[1:EEG.nbchan],3,3,pipeline_visualizations_semiautomated,...
                    0,pipeline_visualizations_semiautomated,[],pipeline_visualizations_semiautomated);
            else
                EEG = pop_eegthresh(EEG,1,[ROI_indices_in_selected_chanlocs]',[reject_min_amp],[reject_max_amp],[EEG.xmin],[EEG.xmax],2,0);
                EEG = pop_jointprob(EEG,1,[ROI_indices_in_selected_chanlocs]',3,3,pipeline_visualizations_semiautomated,...
                    0,pipeline_visualizations_semiautomated,[],pipeline_visualizations_semiautomated);
            end
            EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
            EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0);
            EEG = eeg_checkset(EEG );
            EEG = pop_saveset(EEG, 'filename',strrep(FileNames{current_file}, src_file_ext,'_segments_postreject.set'),'filepath',savepath{3});
        end
        
        %% interpolate the channels that were flagged as bad earlier:
        EEG = pop_interp(EEG, full_selected_channels, 'spherical');
        EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int';
        EEG = eeg_checkset(EEG );
        
        %% re-reference the data: average reference used here
        
        if average_rereference == 1;
            EEG = pop_reref(EEG, []);
            EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int_avgreref';
            EEG = eeg_checkset(EEG);
        else
            [~,ref_chan_indices_in_full_selected_chanlocs] = intersect({full_selected_channels.labels},NO_AVERAGE_REREF_channel_subset);
            EEG = pop_reref(EEG, ref_chan_indices_in_full_selected_chanlocs);
            EEG.setname='wavcleanedEEG_ICA_MARA_rej_chan_int_chansubsetreref';
            EEG = eeg_checkset(EEG);
        end
        
        %% store outputs and report metrics
        Number_Channels_User_Selected(current_file)=size(chan_IDs,2);
        Number_ICs_Rejected(current_file)=length(artifact_ICs);
        Number_Good_Channels_Selected(current_file)=size(selected_channel_locations,2);
        Percent_Good_Channels_Selected(current_file)=Number_Good_Channels_Selected(current_file)/Number_Channels_User_Selected(current_file)* 100;
        Percent_ICs_Rejected(current_file)=Number_ICs_Rejected(current_file)/Number_Good_Channels_Selected(current_file)* 100;
        Percent_Variance_Kept_of_Post_Waveleted_Data(current_file)=varianceWav;
        if isempty(bad_channels_removed)
            Interpolated_Channel_IDs{current_file} = 'none';
        else
            Interpolated_Channel_IDs{current_file}=[sprintf('%s ',bad_channels_removed{1:end-1}),bad_channels_removed{end}];
        end
        Median_Artifact_Probability_of_Kept_ICs(current_file)=median_artif_prob_good_ICs;
        Mean_Artifact_Probability_of_Kept_ICs(current_file)=mean_artif_prob_good_ICs;
        Range_Artifact_Probability_of_Kept_ICs(current_file)=range_artif_prob_good_ICs;
        Min_Artifact_Probability_of_Kept_ICs(current_file)=min_artif_prob_good_ICs;
        Max_Artifact_Probability_of_Kept_ICs(current_file)=max_artif_prob_good_ICs;
        Number_Segments_Post_Segment_Rejection(current_file)=EEG.trials;
        
        %% save preprocessed dataset with subject ID as either txt file (user specified) or eeglab .set file
        switch save_as_format
            case 1 % txt file
                pop_export(EEG,strrep(FileNames{current_file}, src_file_ext,'_processed.txt'),'transpose','on','precision',8);
            case 2 % .mat file
                save([savepath{4} filesep strrep(FileNames{current_file}, src_file_ext,'_processed.mat')], 'EEG');
            case 3 % .set file
                EEG = pop_saveset(EEG, 'filename',strrep(FileNames{current_file}, src_file_ext,'_processed.set'),'filepath',savepath{4});
        end
        
        %% generate power spectrum and topoplot visualization if user requested:
        %plot the spectrum across channels to evaluate pipeline performance
        if pipeline_visualizations_semiautomated == 0
            figure; pop_spectopo(EEG, 1, [], 'EEG' , 'freq', [[freq_to_plot]], 'freqrange',[[vis_freq_min] [vis_freq_max]],'electrodes','off');
            saveas (gcf,[savepath{5} filesep strrep(FileNames{current_file}, src_file_ext,'_processedspectrum.jpg')]);
        end
    end
end

%% generate output table in the "preprocessed" subfolder listing the subject file name and relevant variables for assesssing how good/bad that datafile was and how well the pipeline worked
outputtable=table(FileNames',File_Length_In_Secs',Number_Channels_User_Selected',Number_Segments_Post_Segment_Rejection',...
    Number_Good_Channels_Selected', Percent_Good_Channels_Selected', Interpolated_Channel_IDs',Number_ICs_Rejected',...
    Percent_ICs_Rejected', Percent_Variance_Kept_of_Post_Waveleted_Data',Median_Artifact_Probability_of_Kept_ICs',...
    Mean_Artifact_Probability_of_Kept_ICs',Range_Artifact_Probability_of_Kept_ICs',Min_Artifact_Probability_of_Kept_ICs',...
    Max_Artifact_Probability_of_Kept_ICs');
outputtable.Properties.VariableNames ={'FileNames','File_Length_In_Secs','Number_Channels_User_Selected','Number_Segments_Post_Segment_Rejection',...
    'Number_Good_Channels_Selected', 'Percent_Good_Channels_Selected', 'Interpolated_Channel_IDs','Number_ICs_Rejected',...
    'Percent_ICs_Rejected', 'Percent_Variance_Kept_of_Post_Waveleted_Data','Median_Artifact_Probability_of_Kept_ICs',...
    'Mean_Artifact_Probability_of_Kept_ICs','Range_Artifact_Probability_of_Kept_ICs','Min_Artifact_Probability_of_Kept_ICs',...
    'Max_Artifact_Probability_of_Kept_ICs'};

rmpath(genpath(cleanline_path));
writetable(outputtable, [savepath{4} filesep allsubj{s} 'HAPPE_all_block_output_table ',datestr(now,'dd-mm-yyyy'),'.csv']);

