% Continuous Dots paradigm 
% Subject remains fixated at all times, and clicks left button with left hand if he/she detects coherent motion to the
% left, or clicks right button with right hand if he/she detects motion to the right.
% Dots are alone on the screen - no saccade targets!

reset(RandStream.getDefaultStream,sum(100*clock)); % to avoid randomisation problem

% ***************************************************** BASIC SET - UP 
clear
close all
clc
commandwindow

SITE ='M'; %M=monash
if SITE=='M'
    TheUsualParamsCRT_Monash % 
elseif SITE == 'T'%Trinity College 
    TheUsualParamsCRT_TCD % this script defines some useful parameters of the monitor, trigger port, etc common to all experiments
end

load parDC2  

%%%%%%%%%% IMPORTANT SETTINGS (note many of the time values in ms won't work out exactly - depends on refresh rate & flicker - it'll be the closest it can be, should be checked in EEG)
par.videoFrate = 120;   % Monitor refresh rate
flickerfreq = [20 20];
% If these are the flicker frequencies used the last time, then reverse them left-right
% if all(par.FlickFdots==flickerfreq) | all(par.FlickFdots==fliplr(flickerfreq))
%     par.FlickFdots = fliplr(par.FlickFdots);      % Flicker frequency in Hz for the dots - one for every patch
% else
    par.FlickFdots = flickerfreq;    % [patch-on-left patch-on-right]
% end

par.numPatches = 2;


par.BGcolor=0;
par.cohLevels = [50]; %50 coherence level for DAT1 study (%)
par.secs_btw_targs = [3.06 5.17 7.29]; % in sec
par.targetDur = 1.88;    % coherent motion duration in seconds
par.dotspeed = 6;       % in degrees per second
par.dotsize = 6;   % in pixels
par.numdots = 150;  % this will be the number in a square; the ones outside the circle will be taken away
par.dotpatchsize = 8;   % degrees of visual angle - diameter of dot patch
% par.patchloc = [0 0];
par.patchloc = [-10 -4; 10 -4;]; % patch location coordinates [x y] in degrees relative to center of screen
par.motionDir = [270];    % in degrees relative to positive x-axis (0=rightward), GL: 0 = right, 180 = left, 270=downward

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.numtargets = 18;    %18 trials per block. 
                        %The number of trials per block must be a factor of: (numPatches x cohLevels x secs_btw_targs (ITI) x motionDir). 
                        %So with 2 patches, 1 cohherence level, 3 ITIs, and 1 motion direction -
                        %2*1*3* = 6. So par.numtargets must be a factor of six. 
                        %Thus 18 is a good number of trials in each block
% par.numPatches*length(par.cohLevels)*length(par.secs_btw_targs)*length(par.motionDir); % numPatches x cohLevels x secs_btw_targs x motionDir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par.rgbGRAY = [221 221 221];

% par.EEGtrigCodes = [25 35 50 70 ; 26 36 51 71];
par.MinTimeBtwClicks = 0.5;     % In sec. When polling for mouse clicks, don't want to be counting the same click any more than once.
par.leadintime = 1000;

dlg_title = 'Dot Motion Task';
while 1
    prompt = {'Enter SUBJECT/RUN/TASK IDENTIFIER:','EEG? (1=yes, 0=no)','Use Eyetracker? (1=yes, 0=no)',...
        'Response hand: (1=left, 2=right)'};
    def = {par.runID,num2str(par.recordEEG),num2str(par.useEL),num2str(par.responsehand)};
%     def = {'g_test','0','0','1','2'}; % just testing on laptop
    answer = inputdlg(prompt,dlg_title,1,def);
    par.runID = answer{1};
    par.recordEEG = str2num(answer{2});
    par.useEL = str2num(answer{3});
    par.responsehand = str2num(answer{4});
    if exist([par.runID '.mat'],'file'), 
        dlg_title = [par.runID '.mat EXISTS ALREADY - CHOOSE ANOTHER, OR DELETE THAT ONE IF IT IS RUBBISH']
    else
        break;
    end
end

% Set up for triggers
% Set up for triggers
if par.recordEEG
    if SITE=='M'
        port = hex2dec('2568'); %DN Change this port address for the Monash computer
        lptwrite(port,0);
    elseif SITE=='T'
        % Parallel Port for triggers - set to zero to start
        % port = hex2dec('2010');
        % port = hex2dec('2568');
        % lptwrite(port,0);
        ioObj=io64;
        status=io64(ioObj);
        port = hex2dec('EC00');
        io64(ioObj,port,0);
    end
end

if par.useEL, ELCalibrateDialog, end

window = Screen('OpenWindow', whichScreen, par.BGcolor);

if abs(hz-par.videoFrate)>1
    error(['The monitor is NOT SET to the desired frame rate of ' num2str(par.videoFrate) ' Hz. Change it.'])
end

%  ************************************************* CODES AND TRIAL SEQUENCE
% trigger codes - can only use these 15: [1 4 5 8 9 12 13 16 17 20 21 24 25 28 29]  (don't ask!)
par.CD_RESP  = 1;
par.CD_DOTS_ON = 5;
par.CD_COHMOTION_ON = 9;
par.CD_BUTTONS = [12 13];   % left and right mouse

%%%%%  coherence, motion direction waveforms
nMD = length(par.motionDir);
nDP = size(par.patchloc,1);     
par.FlickTdots = round(par.videoFrate./par.FlickFdots);      % In number of video frames
nrefON = floor(par.FlickTdots/2);                       % number of refreshes where the dots are ON
nrefOFF = par.FlickTdots-nrefON;                      % number of refreshes where the dots are OFF

numcond = 0;
clear conditiondescrip coh trigger triggerFr correctDir patchWithMo
for c=1:length(par.cohLevels)
    for o=1:length(par.secs_btw_targs)
        for p=1:par.numPatches
            numFr1 = round((par.targetDur+par.secs_btw_targs(o)).*par.FlickFdots(p));
            cohMotionOnFr1 = round(par.secs_btw_targs(o)*par.FlickFdots(p))+1;
            for m=1:nMD
                numcond = numcond+1;
                coh{numcond}(1:nMD,1:numFr1) = 0;
                coh{numcond}(m,cohMotionOnFr1:end) = par.cohLevels(c);
                if SITE=='M'|SITE=='E'
                    trigger{numcond} = [par.CD_DOTS_ON 100+numcond];
                elseif SITE=='C'
                    trigger{numcond} = [par.CD_DOTS_ON par.CD_COHMOTION_ON];
                end
                triggerFr{numcond} = [1 cohMotionOnFr1];
                
                correctDir(numcond) = m;
                patchWithMo(numcond) = p;
                conditiondescrip{numcond} = ['Trigger ' num2str(numcond) ': coherence ' num2str(par.cohLevels(c)) ', motion dir ' num2str(par.motionDir(m)) ', ITI ' num2str(par.secs_btw_targs(o)) ', patch ' num2str(p)];
%                 if par.switchDir 
%                     numcond = numcond+1;
%                     coh{numcond}(1:nMD,1:numFr1) = 0;
%                     coh{numcond}(3-m,cohMotionOnFr1:par.switchDirFr(o)-1) = par.cohLevels(c);
%                     % transition period
%                     si = linspace(0,par.cohLevels(c),par.numswitchFr+2);
%                     % ramp down the wrong direction:
%                     coh{numcond}(3-m,par.switchDirFr(o)+[0:par.numswitchFr]) = fliplr(si(1:end-1));
%                     % 5 ramp up the right direction:
%                     coh{numcond}(m,par.switchDirFr(o)+[0:par.numswitchFr]) = si(2:end);
%                     coh{numcond}(m,par.switchDirFr(o)+par.numswitchFr+1:end) = par.cohLevels(c);
%                     %                 modir{numcond}(1:par.numFr(o)) = par.motionDir(3-m)*pi/180;
%                     %                 modir{numcond}(par.switchDirFr(o)+length(find(si(2:end-1)<0,1))-1:end) = par.motionDir(m)*pi/180;
%                     trigger{numcond} = trigger{numcond-1};  % triggers are the same as for non switch.
%                     triggerFr{numcond} = triggerFr{numcond-1};
%                     correctDir(numcond) = m;
%                     patchWithMo(numcond) = p;
%                     conditiondescrip{numcond} = ['Trigger ' num2str(numcond) ': coherence ' num2str(par.cohLevels(c)) ', patch ' num2str(p) ', motion dir ' num2str(par.motionDir(m)) ', coh motion onset ' num2str(par.cohMotionOnset(o)) ' SWITCH']
%                 end
            end
        end
    end
end

% No colored cues implemented for continuous dots...
for n=1:numcond
    dotcolor{n} = par.rgbGRAY'*ones(1,size(coh{n},2));
    numperminblock(n) = 1;
end

par.conditiondescrip = conditiondescrip;

disp(['Number of conditions: ' num2str(numcond)])
% Trial condition randomization, new GL:
minblock = [];
for n=1:numcond % -1 because of catch trials, could set to -par.catchTrials
    % minblock = 1:24, if numperminblock == 2, minblock =
    % 1,1,2,2,3,3,...24,24
    minblock = [minblock ones(1,numperminblock(n))*n];
end

% % left, right
left_target_logical = [1,3,5];
right_target_logical = [2,4,6];
temp = [];
number_of_minblocks = ceil(par.numtargets/size(minblock,2));
for blocker = 1:number_of_minblocks
    minblock_rand = minblock(randperm(size(minblock,2)));
%     this_is_crap = 1;
%     while this_is_crap
%         minblock_rand = minblock(randperm(size(minblock,2)));
%         this_is_crap = 0;
%         for i = 1:length(minblock_rand)-4
%             if all(ismember(minblock_rand(i:i+4),left_target_logical)) | all(ismember(minblock_rand(i:i+4),right_target_logical))
%                 this_is_crap = 1;
%             end
%         end
%     end            
    temp = [temp,minblock_rand];
end

% if numtargets is not evenly divisible by the minimum block length, shave off some trials from the end
temp(par.numtargets+1:end)=[];
trialCond = temp;
save(['Quicksave_',par.runID],'trialCond')

if par.useEL
    ELsetupCalib
end
% keyboard
% Screen('FillRect',window, par.BGcolor); % screen blank

window = Screen('OpenWindow', whichScreen, par.BGcolor);

% *********************************************************************************** START TASK
% Instructions:
leftmargin = 0.1;

Screen('DrawText', window, 'Keep your eyes on the middle dot at all times.', leftmargin*scres(1), 0.1*scres(2), 255);
if par.responsehand==1
Screen('DrawText', window, 'Click LEFT mouse button with LEFT hand as FAST as you can', leftmargin*scres(1), 0.2*scres(2), 255);
else 
Screen('DrawText', window, 'Click RIGHT mouse button with RIGHT hand as FAST as you can', leftmargin*scres(1), 0.2*scres(2), 255);
end
Screen('DrawText', window, 'as soon as you are sure you see downward motion.', leftmargin*scres(1), 0.3*scres(2), 255);
Screen('DrawText', window, 'The motion could be in ANY patch.', leftmargin*scres(1), 0.5*scres(2), 255);
Screen('DrawText', window, 'Press any button to begin task.', leftmargin*scres(1), 0.9*scres(2), 255);
Screen('Flip', window); 
HideCursor;
% Things that we'll save on a trial by trial basis
clear PTBtrigT PTBtrig ClickT Click RespLR
RespT=[];
nPTBtrig=0;
numResp=1;

% Waits for the user to press a button.
[clicks,x,y,whichButton] = GetClicks(whichScreen,0);
if par.recordEEG, sendtrigger(par.CD_RESP,port,SITE,0), end
if par.useEL     
    Eyelink('Message', ['TASK_START']);
end
ClickT(1) = GetSecs;
Click(1)=whichButton(1);    % The first response will be the one that sets the task going, after subject reads instructions

%%%%%%%%%%%%%%%%%%%% START TRIALS

% initial lead-in:
% Screen('FillRect',window, 255, fixRect);
% Screen('Flip', window);

Screen('DrawText', window, 'Loading...', 0.35*scres(1), 0.5*scres(2), 255);
Screen('Flip', window);
WaitSecs(par.leadintime/1000);

clear PT
% First make ALL dot stimuli and store:
dots = cell(2,par.numtargets);   % This will contain dot locations relative to center of screen, in DEGREES
for n=1:par.numtargets
    % First make incoherent motion for all patches:
    pwm = patchWithMo(trialCond(n));
    for p=1:par.numPatches
        dots{p,n}=[];
        numFr = round(size(coh{trialCond(n)},2)*par.FlickTdots(pwm)/par.FlickTdots(p));
        % First generate dots at random locations on each frame
        for i=1:numFr
            for d=1:par.numdots
                dots{p,n}(d,:,i) = [(rand-0.5)*par.dotpatchsize (rand-0.5)*par.dotpatchsize];
            end
        end
        % if this is the patch with dots, make coherent motion:
        if p==pwm
            % then add the coherence by selecting dots to move in certain direction relative to
            % previous frame. A different random set is selected each frame.
            for i=2:numFr
                r = randperm(par.numdots);
                for m=1:nMD
                    ncd = round(par.numdots*coh{trialCond(n)}(m,i)/100);
                    randsel = r(1:ncd);
                    % for the selected dots, move them in a particular direction
                    dots{pwm,n}(randsel,1,i) = dots{pwm,n}(randsel,1,i-1)+cos(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots(pwm);         % x-coordinate
                    dots{pwm,n}(randsel,2,i) = dots{pwm,n}(randsel,2,i-1)-sin(par.motionDir(m)*pi/180)*par.dotspeed/par.FlickFdots(pwm);         % y-coordinate
                    r(1:ncd)=[];
                end
                % if it's gone off to the left, wrap it around to the far right
                dots{pwm,n}(find(dots{pwm,n}(:,1,i)<par.dotpatchsize/2),1,i) = dots{pwm,n}(find(dots{pwm,n}(:,1,i)<par.dotpatchsize/2),1,i)+par.dotpatchsize;
                % if it's gone off to the right, wrap it around to the far left
                dots{pwm,n}(find(dots{pwm,n}(:,1,i)>par.dotpatchsize/2),1,i) = dots{pwm,n}(find(dots{pwm,n}(:,1,i)>par.dotpatchsize/2),1,i)-par.dotpatchsize;
                % if it's gone off to the left, wrap it around to the far right
                dots{pwm,n}(find(dots{pwm,n}(:,2,i)<par.dotpatchsize/2),2,i) = dots{pwm,n}(find(dots{pwm,n}(:,2,i)<par.dotpatchsize/2),2,i)+par.dotpatchsize;
                % if it's gone off to the right, wrap it around to the far left
                dots{pwm,n}(find(dots{pwm,n}(:,2,i)>par.dotpatchsize/2),2,i) = dots{pwm,n}(find(dots{pwm,n}(:,2,i)>par.dotpatchsize/2),2,i)-par.dotpatchsize;
            end
        end
        % Finally, go through the dots and get rid of the dots falling outside the
        % circle - put them off the screen.
        for i=1:numFr
            for d=1:par.numdots
                if sqrt(sum(dots{p,n}(d,:,i).^2)) > par.dotpatchsize/2
                    dots{p,n}(d,:,i) = 2*center/deg2px + 0.01;
                end
            end
        end
        PT1 = [];
        for i=1:numFr
            PT1 = [PT1 ; i*ones(nrefON(p),1);zeros(nrefOFF(p),1)];
        end
        PT{n}(:,p) = PT1;
    end   

end

% initial lead-in:
Screen('FillRect',window, 255, fixRect);
Screen('Flip', window);
WaitSecs(3);

% START STIMULATION
portUP=0; lastTTL=0; ButtonDown=0;
for n=1:par.numtargets
    pwm = patchWithMo(trialCond(n));
    % DOT MOTION
    trigs_sent = 0;
    for i=1:size(PT{n},1)
        if par.recordEEG, if SITE=='M'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
        for p=1:par.numPatches
            if PT{n}(i,p)
                Screen('DrawDots', window, dots{p,n}(:,:,PT{n}(i,p))'*deg2px, par.dotsize, dotcolor{trialCond(n)}(:,min([PT{n}(i,p) size(dotcolor{trialCond(n)},2)])), round(center+par.patchloc(p,:).*[1 -1]*deg2px));
            end
        end
        Screen('FillRect',window, 255, fixRect);
            
        trg = find(triggerFr{trialCond(n)}==PT{n}(i,pwm));
        if ~isempty(trg)
            if trg>trigs_sent
                if par.recordEEG, sendtrigger(trigger{trialCond(n)}(trg),port,SITE,1); portUP=1; lastTTL=GetSecs; end
                if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) '_' num2str(trigger{trialCond(n)}(trg))]); end
                nPTBtrig = nPTBtrig+1;
                [VBLTimestamp PTBtrigT(nPTBtrig)] = Screen('Flip', window);
                PTBtrig(nPTBtrig) = trigger{trialCond(n)}(trg);
                trigs_sent = trigs_sent+1;
            else
                Screen('Flip', window);
            end
            checkButton
        else
            Screen('Flip', window);
            checkButton
        end
    end
end
% And then finish out with half of the first trial (only incoherent bit),
% so it doesn't end on a target
n=1;
pwm = patchWithMo(trialCond(n));
% DOT MOTION
trigs_sent = 0;
for i=1:round(size(PT{n},1))
    if par.recordEEG, if SITE=='M'|SITE=='E', if portUP & GetSecs-lastTTL>0.01, lptwrite(port,0); portUP=0; end, end, end
    for p=1:par.numPatches
        if PT{n}(i,p)
            Screen('DrawDots', window, dots{p,n}(:,:,PT{n}(i,p))'*deg2px, par.dotsize, dotcolor{trialCond(n)}(:,min([PT{n}(i,p) size(dotcolor{trialCond(n)},2)])), round(center+par.patchloc(p,:).*[1 -1]*deg2px));
        end
    end
    Screen('FillRect',window, 255, fixRect);

    trg = find(triggerFr{trialCond(n)}==PT{n}(i,pwm));
    if ~isempty(trg)
        if trg>trigs_sent
            if par.recordEEG, sendtrigger(trigger{trialCond(n)}(trg),port,SITE,1); portUP=1; lastTTL=GetSecs; end
            if par.useEL, Eyelink('Message', ['TRIAL' num2str(n) '_' num2str(trigger{trialCond(n)}(trg))]); end
            nPTBtrig = nPTBtrig+1;
            [VBLTimestamp PTBtrigT(nPTBtrig)] = Screen('Flip', window);
            PTBtrig(nPTBtrig) = trigger{trialCond(n)}(trg);
            trigs_sent = trigs_sent+1;
        else
            Screen('Flip', window);
        end
        checkButton
    else
        Screen('Flip', window);
        checkButton
    end
end
if par.useEL 
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',[answer{1},'.edf']);
end

% FEEDBACK %DN:added Feedback here, took code from Redmond's 'reverse'scrip
response_deadline = par.targetDur+0.5;  % in seconds, OFFLINE response deadline for counting performance
if SITE=='M'
    cohmo_trigs = find(PTBtrig>100 & PTBtrig<=100+numcond);
elseif SITE=='C'
    cohmo_trigs = find(PTBtrig==par.CD_COHMOTION_ON);
end
if length(cohmo_trigs)-1 ~=length(trialCond), error('trigger number mismatch!!'); end %DN: Made length(cohmo_trigs)-1 instead of justlength(cohmo_trigs) to get rid of that last half trial maybe delete this change?
clear RTs acc
for n=1:(length(cohmo_trigs)-1) %DN: Made length(cohmo_trigs)-1 instead of justlength(cohmo_trigs) to get rid of that last trial maybe delete this change?
    RTs(n) = nan;
    acc(n) = 0;
    stimtime = PTBtrigT(cohmo_trigs(n));
    nextresp = find(RespT>stimtime & RespT<stimtime+response_deadline,1);
    if ~isempty(nextresp)
        RTs(n) = RespT(nextresp) - stimtime;
        if RespLR(nextresp)==par.responsehand %correctDir(trialCond(n)) %This checks that participant was using the correct hand; 
            acc(n)=1;
        end
    end
end
figure; hist(RTs*1000,[0:100:3000])
title(['MEDIAN RT: ' num2str(median(RTs)*1000)])

disp(['Accuracy: ' num2str(round(mean(acc)*100)) '%'])
disp([num2str(length(find(isnan(RTs)))) ' Misses'])
txtstart = 0.1*scres(1);
Screen('DrawText', window, ['On this block you got ' num2str(round(mean(acc)*100)) '% correct.'], txtstart, 0.25*scres(2), 255);
Screen('DrawText', window, ['You entirely MISSED ' num2str(length(find(isnan(RTs)))) ' trials out of the ' num2str(length(RTs)) '.'], txtstart, 0.35*scres(2), 255);
Screen('DrawText', window, ['Click to Exit'], txtstart, 0.55*scres(2), 255);
Screen('Flip', window); 
GetClicks(whichScreen,0);

%DN: just did this to calculate mean left and right RT and accuracy:
%     LeftTarget_logical = (trialCond == 1|trialCond ==2|trialCond ==5|trialCond ==6|trialCond ==9 ...
%         |trialCond ==10|trialCond ==13|trialCond ==14|trialCond ==17|trialCond ==18|trialCond ==21|trialCond ==22);
%     RightTarget_logical = (trialCond == 3|trialCond ==4|trialCond ==7|trialCond ==8|trialCond ==11 ...
%         |trialCond ==12|trialCond ==15|trialCond ==16|trialCond ==19|trialCond ==20|trialCond ==23|trialCond ==24);
%     LeftTargetRTs = RTs(LeftTarget_logical);
%     RightTargetRTs = RTs(RightTarget_logical);
%     MeanLeftTargetRT=nanmean(LeftTargetRTs)
%     MeanRightTargetRT=nanmean(RightTargetRTs)
%     %%%%%%%%Accuracy
%     LeftTargetAcc = acc(LeftTarget_logical);
%     RightTargetAcc = acc(RightTarget_logical);
%     LeftTargetAccPercent = mean(LeftTargetAcc)*100;  
%     RightTargetAccPercent = mean(RightTargetAcc)*100;

%DN: now saves mean left and right RT and accuracy
%     save([par.runID],'ClickT','Click','nPTBtrig','PTBtrigT','PTBtrig','RespT','RespLR','trialCond','coh','dotcolor','trigger','triggerFr','correctDir','par', 'conditiondescrip','RTs','acc','MeanLeftTargetRT','MeanRightTargetRT','LeftTargetAccPercent','RightTargetAccPercent')
responsehand=par.responsehand;
save([par.runID],'ClickT','Click','nPTBtrig','PTBtrigT','PTBtrig','RespT','RespLR','trialCond','coh','dotcolor','trigger','triggerFr','correctDir','par', 'conditiondescrip','RTs','acc', 'responsehand')
save parDC2 par

cleanup