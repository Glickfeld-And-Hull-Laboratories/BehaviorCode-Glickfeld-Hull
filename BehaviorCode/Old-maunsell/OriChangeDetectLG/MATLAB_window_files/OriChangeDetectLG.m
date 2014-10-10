function [retval] = HoldAndDetectConstant(data_struct, input)

% This function processes event data and saves variables in input.
%
% Then it calls plot functions in a
% try-catch block so that errors in plot functions don't affect
% variable saving.
%  MH 100115
% $Id$
if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end
beep off;  % otherwise played through speakers to animal!
format compact; 

ds = data_struct;

addpath('~/Repositories/MWorksMatlabToolbox');
addpath('~/Repositories/tools-mh');

% debug
%dumpEvents(ds); 
%length(ds.events)
%save('/tmp/test.mat', 'data_struct', 'input');

% First call after pressing "play" in client
if nargin == 1 || ~isfield(input, 'trialSinceReset')
    % init input
    disp('First trial, initializing input');
    input.trialSinceReset = 1;
    input.startDateVec = datevec(now);
    input.reactTimesMs = {};
    input.holdTimesMs = {};
    input.reqHoldTimeMs = {};
    input.tooFastTimeMs = [];
    input.trialOutcomeCell = {};
    input.holdStartsMs = {};
    input.juiceTimesMsCell = {};
    input.constChangedOnTrial = {};
    input.laserPowerNum = {};
    input.laserPowerMw = {};
    input.savedEvents = {};
    input.eventCodecs = {};
    input.gratingDurationMs = {};
    input.laserPeakMaxMw = {};
    input.laserRampLengthMs = {};
    input.tItiWaitTimeMs = {};
    input.tGratingDurationMs = {};
    input.tLaserPeakMaxMw = {};
    input.tLaserRampLengthMs = {};
    input.tLaserDoLinearRamp = {};
    input.tLaserDoPulseTrain = {};
    input.tBlock2TrialNumber = {};
    input.tGotLongBonusThisTrial = {};
    input.tTrialLaserPowerMw = {};
    input.tTrialLaserOnTimeMs = {};
    input.tTrialLaserOffTimeMs = {};
    input.tGratingElevationDeg = {};
    input.tGratingAzimuthDeg = {};
    input.tGratingHeightDeg = {};
    input.tGratingWidthDeg = {};
    input.tGratingSpatialFreqCPD = {};
    input.tGratingSpeedDPS = {};
    input.tGratingDurationMs = {};
    input.tTargetGratingAlphaMultiplier = {};
    input.tBaseGratingAlphaMultiplier = {};
    input.tGratingBaseDirection = {};
    input.tDirectionStepDeg = {};
    input.tGratingDirectionDeg = {};
    input.tGratingDirectionStepsPerOctave = {};
else
  input.trialSinceReset = input.trialSinceReset+1;
  assert(all(input.tooFastTimeMs >= 1));
end
thisTrialN = input.trialSinceReset;

%% quick path check 
if ~exist('mwGetEventValue')
  error('Missing mwGetEventValue - is path correct?');
end

events = ds.events;

%% save codec
if thisTrialN == 1
  input.eventCodecs{1} = ds.event_codec;
else
  lastCodecN = find(cellfun(@(x) ~isempty(x), input.eventCodecs), ...
                    1, 'last');
  if ~isequalwithequalnans(input.eventCodecs{lastCodecN}, ds.event_codec)
    error('bug: Codec seems to have changed - stopping!');
    %input.eventCodecs{thisTrialN} = ds.event_codec;
  end
end


%% process #announceStimulus events

asCodeN = codec_tag2code(input.eventCodecs{1}, '#announceStimulus');
codeV = cat(2, ds.events.event_code);
timeV = cat(2, ds.events.time_us);
asIx = codeV == asCodeN;
asTimeV = timeV(asIx);
%ds.events(asIx).data;

allData = {ds.events(asIx).data};
asActionV = cellfun(@(x) x.action, allData, 'UniformOutput', false);
asNameV = cellfun(@(x) x.name, allData, 'UniformOutput', false);
assert(all(strcmp(asActionV, 'draw')), ...
       'Found an action that was not ''draw'' - fix!');
unqNames = unique(asNameV);
for iN = 1:length(unqNames)
  tName = unqNames{iN};
  tNameIx = strcmp(asNameV, tName);
  tTimes = asTimeV(tNameIx);
  input.announceStimulusTimes(1).(tName){thisTrialN} = tTimes;
end
% remove  from savedEvents
events = events(~asIx);

% display some data about it
% $$$ at0 = input.announceStimulusTimes;
% $$$ wt0 = whos('at0');
% $$$ at1 = input.savedEvents;
% $$$ wt1 = whos('at1');
% $$$ fprintf(1, '--------------- savedEvents %.2f MB, asTimes %.2f MB\n', ...
% $$$         wt1.bytes/1e6, wt0.bytes/1e6);




%% save event stream
input.savedEvents{thisTrialN} = events;

%% are constants here?  if so, extract
nSync = length(mwGetEventValue(events, ds.event_codec, ...
                             'sync', 'all', 'ignoreMissing'));
constsStruct = [];
allEvents = events;
if nSync >= 4 % is Justin Timberlake missing?
  % constants arrived.  Split out
  disp('Found dump of constants, splitting out from events');
  
  syncCode = codec_tag2code(ds.event_codec, 'sync');
  codeList = [events.event_code];
  syncNs = find(codeList == syncCode);

  % check a lot of assumptions.
  %   What should happen is sync=1, constants dumped included sync
  %   = 1, first trial starts, is interrupted by pressing stop,
  %   restart sends new codes with sync = 1, then sync = 0 stops
  %   trial.  Don't ask me why this is, pressing it once doesn't
  %   send the constants.
  errFlag = false;
  if ~ all(diff([events.time_us]) >= 0)
    disp('increasing in time?');     errFlag = true;
  elseif ~(syncNs(1) == 1)
    disp('sync structure'); errFlag = true;
  elseif ~(syncNs(end) == length(events))
    disp('last code'); errFlag = true;
  elseif ~(events(syncNs(1)).data == 1), 
    disp('err1'); errFlag = true; 
  elseif ~(events(syncNs(2)).data == 1), 
    disp('err2'); errFlag = true;
  elseif ~(events(syncNs(end-1)).data == 1), 
    disp('err3'); errFlag = true;
  elseif ~(events(syncNs(end)).data == 0), 
    disp('err4'); errFlag = true;
  end

  if errFlag == true
    disp('** error in sync code pattern: check and fix code');
    disp(syncNs);
    disp(nSync);
    retval = []
    return
  end
  
  if nSync > 4
    disp('Too many sync codes:');
    disp(nSync);
    % just display a message
    %   errFlag = true;
  end

  % if 6 codes, pull out consts from first block and leave last
  % (Somewhat dangerous)
  constsStruct = events(1:syncNs(3));
  trialStruct = events(syncNs(end-1):end);
  %events = trialStruct;
  
  disp(sprintf('** Constants arriving in codec this trial! (tr %d)', ...
               input.trialSinceReset));
else
  trialStruct = events;
end


%% make sure we're getting constants on first trial
if thisTrialN == 1
  errFlag = false;
  if isempty(constsStruct)
    errFlag = true;
  else
    tSubjectNumTemp = mwGetEventValue(constsStruct, ds.event_codec, 'subjectNum', ...
                                      [], 'ignoremissing');
    if isempty(tSubjectNumTemp)
      errFlag = true;
    end
  end
  if errFlag
    disp('** Error: must press start button twice on first trial!');
    retval = [];
    return
  end
end

%% process constants
input.constList = ...  % all the saved variables
    { 'doExtendItiOnShortPrevTrial', ...
      'doLaserStim', ...
      'doLever', 'doLeverSolenoidAllTrials', ...
      'doLeverSolenoidOnEarly', 'doLeverSolenoidOnMiss', ...
      'doVisualStim', ...
      'doGeomHoldDist', ...
      'geomHoldMeanMs', ...
      'earlyTimeoutMs', 'experimentXmlTrialId', ...
      'fakeMouseMaxPressMs', 'fakeMouseMaxReactMs', ...
      'fixedReqHoldTimeMs', 'forceMinStimStepTo', ...
      'targetGratingContrast', 'baseGratingContrast', ...
      'gratingBaseDirection', 'gratingDirectionStepDeg','gratingDirectionStepsPerOctave', ...
      'gratingAzimuthDeg', 'gratingElevationDeg',...
      'gratingWidthDeg','gratingHeightDeg', ... 
      'gratingSpatialFreqCPD', ...
      'gratingSpeedDPS',  ...
      'gratingDurationMs', ...
      'interRewardIntervalMs', 'itiTimeMs', ...
      'itiExtraRandTimeMs', ...
      'laserDoLinearRamp', 'laserDoPulseTrain', ...
      'laserOffPowerMw', 'laserPeakMaxMw', ...
      'laserPeakStepsPerOctave', 'laserPulseLengthMs', ...
      'laserPulsePeriodMs', 'laserRampLengthMs', ...
      'laserTrainLengthMs', 'laserTransitionRampUpDownMs', ...
      'laserTransitionDoExpRamp', ...
      'laserRampExtraConstantLengthMs', ...
      'laserRampDoExpRamp', ...
      'trialLaserPowerMw', ...
      'trialLaserOnTimeMs', ...
      'trialLaserOffTimeMs', ...
      'maxConsecCorrects', ...
      'nConsecErrorsCauseTimeout', ...
      'consecErrorTimeoutS', ...
      'doLongBonus', ...
      'longBonusExtraHoldTimeMs', ...
      'longBonusExtraRewardUs', ...
      'maxRewardUs', ...
      'minRewardUs', ...
      'missedTimeoutMs', ...
      'postRewardMs', ...  
      'randReqHoldMaxMs', ...
      'reactTimeMs', ...
      'subjectNum', ...
      'trPer80Level1', 'trPer80Level2', 'trPer80Level3', ...
      'trPer80Level4', 'trPer80Level5', 'trPer80Level6', ...
      'trPer80Level7', 'trPer80Level8', ...
      'doBlock2', ...
      'block2DoRampLength', 'block2RampLengthMs2', ...
      'block2RampLengthPowerMaxMw2', ...
      'block2DoRampVTrain', ...
      'block2RvtTrainPowerMaxMw', ...
      'block2RvtTrainStepsPerOctave', ...
      'block2DoGratingAppearance', ...
      'block2TargetGratingContrast', ...
      'block2GratingBaseDirection', 'block2GratingDirectionStepDeg', 'block2GratingDirectionStepsPerOctave', ...
      'block2GratingElevationDeg', ...
      'block2GratingAzimuthDeg', ...
      'block2GratingHeightDeg', ...
      'block2GratingWidthDeg', ...
      'block2GratingSpatialFreqCPD', ...
      'block2GratingSpeedDPS', ...
      'block2GratingDurationMs', ...
      'block2DoRampVTrain', 'block2RvtTrainPowerMaxMw', ...
      'block2DoTrialLaser', 'block2TrialLaserPowerMw', ...
      'block2TrialLaserOnTimeMs', ...
      'block2TrialLaserOffTimeMs', ...
      'block2TrialLaserTargetGratingContrast', 'block2TrialLaserBaseDirection', 'block2TrialLaserDirectionStepDeg', 'block2TrialLaserDirectionStepsPerOctave', ...
      'tooFastTimeMs' };

nConsts = length(input.constList);

% look for changes, all trials  (events)
b={};
fieldList = fieldnames(input);
if input.trialSinceReset > 1
  for iC = 1:nConsts
    tCN = input.constList{iC};
    tV = mwGetEventValue(allEvents, ds.event_codec, tCN, 'last', 'ignoremissing');
    
    if ~isempty(tV)
      if ~any(ismember(fieldList, tCN)) || ~isequalwithequalnans(tV, input.(tCN))
        % different than last trial (or not a field - only happens when adding vars to xml)
        % stuff into array saying what changed
        b(end+1,1:2) = {tCN, tV};
      end
    end
  end
end
input.constChangedOnTrial{input.trialSinceReset} = b;


% get constants from dump, first trial only  (constsStruct)
for iC = 1:nConsts
  tCN = input.constList{iC};
  if input.trialSinceReset==1  % use dump on first trial
    tCV1 = mwGetEventValue(constsStruct, ds.event_codec, tCN, [], ...
                          'ignoremissing');

    if isempty(tCV1)
      error(sprintf('*** Constant %s not found in dump - is it checked in Events box?!!!', ...
                   tCN));
    end
    input.(tCN) = tCV1;
  
  else  % if not trial 1, always copy const values in
    tCV = mwGetEventValue(events, ds.event_codec, tCN, 'last', ...
                          'ignoremissing');
    if ~isempty(tCV)
      input.(tCN) = tCV;
    end
  end
end

% make some into arrays
input.trPer80V = [input.trPer80Level1 input.trPer80Level2 input.trPer80Level3 ...
                  input.trPer80Level4 input.trPer80Level5 input.trPer80Level6 ...
                  input.trPer80Level7 input.trPer80Level8];


% save first trial status
if input.trialSinceReset == 1
  input.firstTrConsts = input;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% process corrects
nS = length(mwGetEventValue(trialStruct, ds.event_codec, ...
                             'success', 'all', 'ignoreMissing'));
nF = length(mwGetEventValue(trialStruct, ds.event_codec, ...
                             'failure', 'all', 'ignoreMissing'));
nI = length(mwGetEventValue(trialStruct, ds.event_codec, ...
                             'ignore', 'all', 'ignoreMissing'));
if sum(nS+nF+nI) > 1, 
  disp(trialStruct); disp(ds.event_codec); disp(nSync);
  disp(nS); disp(nF); disp(nI);
  error('error - 1, check'); 
end
if nS == 1, tStr = 'success'; 
elseif nF == 1, tStr = 'failure';
elseif nI == 1, tStr = 'ignore';
else
  disp('Error!  Missing trial outcome variable this trial');
  tStr = 'error-missing';
end
% save
input.trialOutcomeCell{trN} = tStr;

%% process reaction times for this trial
codes = [trialStruct.event_code];


stimOnUs = mwGetEventTime(trialStruct, ds.event_codec, 'stimulusOn', 1);
totalHoldTimeMs = mwGetEventValue(trialStruct, ds.event_codec, 'tTotalReqHoldTimeMs');
%trialStartUs = mwGetEventTime(trialStruct, ds.event_codec, 'trialStart', 1);
leverDownUs = mwGetEventTime(trialStruct, ds.event_codec, 'leverResult', 1, 1);
leverUpUs = mwGetEventTime(trialStruct, ds.event_codec, 'leverResult', 2, 0);

%reactTimeMs = (leverUpUs - stimOnUs) / 1000;
holdTimeMs = (leverUpUs - leverDownUs) / 1000;
reactTimeMs = (holdTimeMs - totalHoldTimeMs);

% add to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;
input.reqHoldTimeMs{trN} = totalHoldTimeMs;

input.tItiWaitTimeMs{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tItiWaitTimeMs', 'last', 1);

% add to array: block2 settings
input.tGratingDurationMs{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingDurationMs', 'last', 1);
input.tLaserPeakMaxMw{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tLaserPeakMaxMw', 'last', 1);
input.tLaserRampLengthMs{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tLaserRampLengthMs', 'last', 1);
input.tBlock2TrialNumber{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tBlock2TrialNumber', 'last', 1);
input.tLaserDoLinearRamp{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tLaserDoLinearRamp', 'last', 1);
input.tLaserDoPulseTrain{trN}= mwGetEventValue(trialStruct, ds.event_codec, 'tLaserDoPulseTrain', 'last', 1);
input.tTrialLaserPowerMw{trN}= mwGetEventValue(trialStruct, ds.event_codec, 'tTrialLaserPowerMw', 'last', 1);
input.tTrialLaserOnTimeMs{trN}= mwGetEventValue(trialStruct, ds.event_codec, 'tTrialLaserOnTimeMs', 'last', 1);
input.tTrialLaserOffTimeMs{trN}= mwGetEventValue(trialStruct, ds.event_codec, 'tTrialLaserOffTimeMs', 'last', 1);
input.tGratingDirectionDeg{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingDirectionDeg', 'last', 1);
input.tGratingAlphaMultiplier{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tTargetGratingAlphaMultiplier', 'last', 1);
input.tGratingAlphaMultiplier{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tBaseGratingAlphaMultiplier', 'last', 1);
input.tGratingBaseDirection{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingBaseDirection', 'last', 1);
input.tGratingBaseDirection{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingDirectionStepsPerOctave', 'last', 1);
input.tGratingElevationDeg{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingElevationDeg', 'last', 1);
input.tGratingAzimuthDeg{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingAzimuthDeg', 'last', 1);
input.tGratingHeightDeg{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingHeightDeg', 'last', 1);
input.tGratingWidthDeg{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingWidthDeg', 'last', 1);
input.tGratingSpatialFreqCPD{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingSpatialFreqCPD', 'last', 1);
input.tGratingSpeedDPS{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingSpeedDPS', 'last', 1);
input.tGratingDurationMs{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingDurationMs', 'last', 1);

% total reward times
juiceAmtsUs = mwGetEventValue(trialStruct, ds.event_codec, 'juice', 'all');
juiceAmtsMs = double(juiceAmtsUs(juiceAmtsUs~=0))/1000;
input.juiceTimesMsCell{thisTrialN} = juiceAmtsMs;

input.tGotLongBonusThisTrial{trN} = mwGetEventValue(trialStruct, ds.event_codec, 'tGotLongBonusThisTrial', 'last', 1);


% array processing
input.trPer80History{trN} = input.trPer80V;

if input.doVisualStim == 1
  % direction
  tGratingDirectionDeg = mwGetEventValue(trialStruct, ds.event_codec, 'tGratingDirectionDeg', 'all');  
  input.gratingDirectionDeg{thisTrialN} = tGratingDirectionDeg;
end
if input.doLaserStim == 1
  % laser power values  
  tLaserPowerMw = mwGetEventValue(trialStruct, ds.event_codec, ...
                                  'tLaserPowerMw', 'all');
  input.laserPowerMw{thisTrialN} = tLaserPowerMw(end);
  %disp(sprintf('Power was %3.2fmW', input.laserPowerMw{thisTrialN}));
end



%% disp status
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
stimStr = '';
if input.doVisualStim
  stimStr = [stimStr sprintf('direction %gdeg ', chop(input.gratingDirectionDeg{thisTrialN},2))];
end
if input.doLaserStim
  stimStr = [stimStr sprintf('power %gmW ', ...
                    chop(input.laserPowerMw{thisTrialN}, 2))];
end
if ~input.doBlock2
  block2Str = '';
else
  block2Str = sprintf('b2tr %d ', input.tBlock2TrialNumber{thisTrialN});
end
if input.itiExtraRandTimeMs > 0
  itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{thisTrialN}));
else
  itiStr = '';
end
fprintf(1,'Hold %d, req %d, react %d, %s%s %s- %d rew %dms\n', ...
        round(holdTimeMs), round(totalHoldTimeMs), round(reactTimeMs), ...
        itiStr, ...
        stimStr, block2Str, ...
        nJ, round(juiceD));
        
%% save variables for next trial
retval = input;

%% run subfunctions
doProfile = false;
if doProfile
  profile on;
end

tic
try
  input = saveMatlabState(data_struct, input);

  % do the plot
  if ~isfield(input, 'lastPlotElapsedS')
    input.lastPlotElapsedS = 0;
  end
  if input.lastPlotElapsedS > 2
    % skip every other plot to make sure we catch up
    disp(sprintf('Skipping plot, last elapsed was %5.1fs', ...
                 input.lastPlotElapsedS));
        input.lastPlotElapsedS = 0;
  else
    % last one was fast, plot this time
    tic
    input = plotOnlineHist(ds, input);  % this gets the CHANGED event stream w/ first trial consts extracted
    input.lastPlotElapsedS = toc;
    %disp(input.lastPlotElapsedS);
  end
    

  %input = testUpload(data_struct, input);

catch ex
  disp('??? Error in subfunction; still saving variables for next trial')
  printErrorStack(ex);
end
tEl = toc;
%if tEl > 1.3
%  disp(sprintf('          Elapsed: %gsec  ', chop(tEl,2)));
%end

if doProfile
  p = profile('info');
  save('/tmp/profileHNDC.mat', p)
end




%% save variables for next trial
retval = input;


return

