function [retval] = CerebellarStim(data_struct, input)
% Online Matlab processing for cerebellarStim.xml
% Code for Hull Lab - Duke Neurobiology
% Completed August 29th, 2014 by Andrew McKinney
% call: exptSetupBridge, exptProcessBridgeInput, package data, then call subfunctions

if nargin < 2,
    input = [];
    input.trialSinceReset = 0;
    input.saveTime = datestr(now, 'yymmdd-HHMM');
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
oneValEachTrialNames = { ...
    'tTempStim', ...
    'tCounter2', ...
    'ttCounter2', ...
    'tttCounter2', ...
    'tTactileStimulusDurationUs', ...
    'spCounter1', ...
    'spCounter2', ...
    'spCounter3', ...
    'spCounter4', ...
    'spCounter5', ...
    'spCounter6', ...
    'spCounter7', ...
    'spCounter8', ...
    'spCounter9', ...
    'spCounter10', ...
    'tITIWheelCounter', ...
    'tISIWheelCounter', ...
    'tStimWheelCounter', ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tStimTurnedOn'};

if isempty(input),
    input.trialSinceReset = 1;
    input.startDateVec = datevec(now);
    input.savedEvents = {};
    input.eventCodecs = {};
    input.eventCodecs{1} = ds.event_codec;
    nOne = length(oneValEachTrialNames);
    for iV = 1:nOne
        input.(oneValEachTrialNames{iV}) = {};
    end
else
    input.trialSinceReset = input.trialSinceReset+1;
end
trN = input.trialSinceReset;

% Consts that govern tValues
input.constList = { 'subjectNum',...
    'experimentXmlTrialId', ...
    'nScansOn', ...
    'speedTimerIntervalMs', ...
    'nScansOff', ...
    'stopAfterNTrials', ...
    'tactileStimulusDurationUs',...
    'soundDurationMs',...
    'postSoundPauseDurationMs',...
    'frameImagingRateMs',...
    'frameImagingExposureMs',...
    'frameImagingBinning'};

events = ds.events;
eventsTrial = events;


for iV = 1:length(oneValEachTrialNames)
    tName = oneValEachTrialNames{iV};
    input.(tName){trN} = mwGetEventValue(eventsTrial, ds.event_codec, tName, 'last', 1);
end

nSync = length(mwGetEventValue(ds.events, ds.event_codec, ...
                             'sync', 'all', 'ignoreMissing'));

                         
syncCode = codec_tag2code(ds.event_codec, 'sync');
codeList = [ds.events.event_code];
syncNs = find(codeList == syncCode);                        
if nSync >= 4,
  % if 6 codes, pull out consts from first block and leave last
  % (Somewhat dangerous)
  constsStruct = struct;
  constsStruct = ds.events(1:syncNs(3));
  trialStruct = ds.events(syncNs(end-1):end);
  ds.events = trialStruct;
end


nConsts = length(input.constList);
for iC = 1:nConsts
  tCN = input.constList{iC};
  if input.trialSinceReset==1  % use dump on first trial
    tCV1 = mwGetEventValue(constsStruct, ds.event_codec, tCN, [], ...
                          'ignoremissing');
    input.(tCN) = tCV1;
  
  else  % look for changes on all trials
    tCV = mwGetEventValue(ds.events, ds.event_codec, tCN, [], ...
                          'ignoremissing');
    if ~isempty(tCV)
      input.(tCN) = tCV;
    end
  end
end
 
% Look for changes in constants for all trials
b={};
for iC = 1:nConsts
  tCN = input.constList{iC};
  tV = mwGetEventValue(ds.events, ds.event_codec, tCN, [], 'ignoremissing');
 
  if ~isempty(tV), 
    % Constructs array of changed consts and when they changed
    b(end+1,1:2) = {tCN, tV};
  end
end
input.constChangedOnTrial{input.trialSinceReset} = b;


input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');

disp('This is a trial! I wrote code that works!');

%% Run Plotting Function
input.holdStartsMs = input.tThisTrialStartTimeMs;
%input = exptRunSubfunctions(ds, input, { @plotVisStim });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Returns calculated and passed variables for the next trial.
retval = input;
return

