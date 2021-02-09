function runafew_Nicole_after(number)
%clear all; number=1;
subj=1:75;
LPF=[35 ];
HPF=[0.01 1000];%1000 is HAPPE
[a b c] = ndgrid(subj, LPF, HPF);
allcombinations = [a(:) b(:) c(:)];
conditions = allcombinations(number,:);
subj = conditions(1); LPFused=conditions(2);HPFused = conditions(3);
%subj=6
clearvars -except subj LPFused HPFused

cd ..
Start
cd Process_Matlab
close all
clc

warning off
%% Monash Dyslexia participants:
path_temp = '../Data/';
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
                 'NV75CT','OB50CT','OH08T','OK34DT','DP06DT','RD15DT','SH81CT','SS85CT','TB28CT','TG51CT',...
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
%secondSessionKids={} Check that all participants need to go here %
allblocks = {[1:3 5:12],[1:5 7:12],[1:14],[1:12],[1:11],[1:12],[1:12],[1:7 9:12],[1:12],[1:12],...%10
             [1:12],[1:12],[1:5,7:13],[1:12],[1:12],[1:9 11:12],[1:12],[1:12],[2:6],[1:12],...%20
             [1:12],[1:9 11:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:13],[1:12],...%30
             [1:12],[1:3 5:12],[1:12],[1:2 4:13],[1:10 12:17],[1:12],[1:12],[1:10,12],[1:5 7:9 11:20],[1:12],...%40
             [1:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:8 10 12],[1:9 11],[1:12],[1:12],...%50
             [1:12],[1:12],[1:6 8:13],[1:8 10:12],[1:12],[1:12],[1:12],[2:7 9:12],[1:12],[1:12],...%60
             [1:12],[1:12],[1:12],[1:2,5:13],[1:6],[1:2 4:10],[1:4 6:12],[1:10 12:13],[1:12],[1:12],...%70
             [1:9 12],[1:12],[1:12],[1:12],[1:12]};%skip 3}; 
%temp blocks to check for bad channels       
allblocks = {[1:3 5:12],[1:5 7:12],[1:14],[1:12],[1:11],[1:12],[1:12],[1:7 9:12],[1:12],[1:12],...%10
             [1:12],[1:12],[1:5,7:13],[1:12],[1:12],[1:9 11:12],[1:12],[1:12],[2:6],[1:12],...%20
             [1:12],[1:9 11:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:13],[1:12],...%30
             [1:12],[1:3 5:12],[1:12],[1:2 4:13],[1:10 12:17],[1:12],[1:12],[1:10,12],[1:5 7:9 11:20],[1:12],...%40
             [1:12],[1:12],[1:12],[1:12],[1:12],[1:12],[1:8 10 12],[1:9 11],[1:12],[1:12],...%50
             [1:12],[1:12],[1:6 8:13],[1:8 10:12],[1:12],[1:12],[1:12],[2:7 9:12],[1:12],[1:12],...%60
             [1:12],[1:12],[1:12],[1:2,5:13],[1:6],[1:2 4:10],[1:4 6:12],[1:10 12:13],[1:12],[1:12],...%70
             [1:9 12],[1:12],[1:12],[1:12],[1:12]};%skip 3}; 

         
         %block 

allbadchans = {[],[],[28 39],[28 32],[],[],[],[],[],[],...
               [],[],[],[],[17],[],[17,22],[],[],[41],...
               [],[17],[28,29,22,32,58,43],[],[],[],[17,22],[],[28,41],[22,28],...
               [],[10,17],[17],[],[],[],[],[],[],[]...
               [17,22,26,46],[],[],[],[],[],[],[17,22,46,47,62,63],[],[2,36],[16]...
               [],[],[],[46],[1 4],[31],[17,22,41,46],[41,46],[17,22,41,46],[17 51]...
               [17,28,64],[1,9,19,25,48],[],[],[],[],[17,22,60],[17,22,41,46],[],[17,22,33,36]...
               [],[28],[],[],[]};
%allbadchans = {[],[],[36 28 39],[28 32 41],[41 32],[],[],[42,16,29,30,35],[27],[],...
%               [46 7 3 37 41],[],[],[],[17 28],[20],[17,22,12,16,27,31,1,2,33,36,7,45,46,3,41,42,8,40,11,50,28,32],[16 12 23 64,1 ,2,33,36],[],[41 17],...
%               [17 28 32 22],[17 32],[28,29,22,32,58,43],[],[],[12 41 42],[17,22],[1],[28,41],[22,28,61],...
%               [29 63],[10,17],[17],[],[],[],[],[],[28 32],[],...
%               [17,22,26,46],[],[64 60],[],[63 64 31],[],[],[17,22,46,47,62,63],[],[2,36 16],...
%               [],[],[15 45 64 60],[46 49],[1 4 11 17 28 32 50 8],[31],[17,22,41,46],[41,46],[31,16,17,22,41,46,64],[17 51],...
%               [17,28,64,45,41,42],[1,9,19,25,48],[40,7],[46],[64,32],[],[17,22,60],[17,22,41,46],[51,41],[17,22,33,36],...
%               [28,29,30,31,32,64,60,63],[28],[],[],[]};  
           % when considering exclusion criteria: do not use the frontal
           % electrodes
%Monash_Stationary_EEG = {'KJ1001'}; %List participants tested using the Bellgrove Lab's stationary EEG system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_label = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'}; %maximum number of blocks for all participants

duds = [];  
single_participants = [subj]; %can put 1 participant, or multiple in here. E.g. [3:8] to process participants 3 to 8

%%

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    allblocks([duds]) = [];
    allbadchans([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder([single_participants]);
    allsubj = allsubj([single_participants]);
    allblocks = allblocks([single_participants]);
    allbadchans = allbadchans([single_participants]);
end

%% CSD
E = textread('chans65_monash.asc','%s');
M = ExtractMontage('10-5-System_Mastoids_EGI129.csd',E,subj);  % reading in the montage for the CSD toolbox
% MapMontage(M);
[G_monash_stationary,H_monash_stationary] = GetGH(M);

%%

for s=1:length(allsubj)
    disp(allsubj{s})
    blocks = allblocks{s};
    badchans = allbadchans{s};
    clear paths files matfiles ET_files ET_matfiles; k=0;
    for n=1:length(blocks)
        k=k+1;
            paths{k} = [path_temp subject_folder{s} '/'];
            files{k} = [allsubj{s} filesep allsubj{s} block_label{blocks(n)} '_processed.mat'];
            filesorig{k} = [allsubj{s} block_label{blocks(n)} '.vhdr'];
            matfiles{k} = [path_temp subject_folder{s} '/' allsubj{s} block_label{blocks(n)}, '.mat'];
            ET_files{k}=[path_temp 'Samples_and_Events/' allsubj{s} block_label{blocks(n)}, '.asc'];
            ET_matfiles{k} = [path_temp subject_folder{s} '/' allsubj{s} block_label{blocks(n)}, '_ET.mat'];
    end            

    
    if ismember(subject_folder{s},Dyslexic)
        G_CSD = G_monash_stationary;
        H_CSD = H_monash_stationary;
        dyscAMRM = 'Dyslexic';
        1
        %Monash_Stationary_EEG_preprocess_Dyslexic %Preprocessing for  Bellgrove Lab's stationary EEG system
    elseif ismember(subject_folder{s},AM_control)
        G_CSD = G_monash_stationary;
        H_CSD = H_monash_stationary;
        dyscAMRM = 'AM_control';
        2
        %Monash_Stationary_EEG_preprocess_Aged_matched_control %Preprocessing for Bellgrove Lab's portable EEG system
    elseif ismember(subject_folder{s},RM_control)
        G_CSD = G_monash_stationary;
        H_CSD = H_monash_stationary;
        dyscAMRM = 'RM_control';
        3
        %Monash_Stationary_EEG_preprocess_Reading_matched_control %Preprocessing for Bellgrove Lab's portable EEG system        
    else
        keyboard %the participant needs to be listed above in either Monash_Stationary_EEG or Monash_Portable_EEG
    end
    
    Monash_Stationary_EEG_preprocessAfter
end
