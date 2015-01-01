function [retval] = CerebellarStim(data_struct, input)
% Online Matlab processing for cerebellarStim.xml
% Code for Hull Lab - Duke Neurobiology
% Completed August 29th, 2014 by Andrew McKinney
% call: exptSetupBridge, exptProcessBridgeInput, package data, then call subfunctions

if nargin < 2,
    input = [];
end

ds = data_struct; %added Couunter ttCounter and Counter2
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
oneValEachTrialNames = { ...
    'counter', ...
    'tCounter', ...
    'ttCounter', ...
    'counter2', ...
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
    'itiCounter1', ...
    'itiCounter2', ...
    'itiCounter3', ...
    'itiCounter4', ...
    'itiCounter5', ...
    'itiCounter6', ...
    'itiCounter7', ...
    'itiCounter8', ...
    'itiCounter9', ...
    'itiCounter10', ...
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
    input.saveTime = datestr(now, 'yymmdd-HHMM');
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
    'nScansOff', ...
    'stopAfterNTrials', ...
    'tactileStimulusDurationUs',...
    'soundDurationMs',...
    'preSoundPauseNFrames', ...
    'postSoundPauseNFrames', ...
    'frameImagingFrequencyHz',...
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

try
    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

try
  input.counter2TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter2', 'all', [], 1);
  input.counter2Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter2', 'all', 1) ;
catch
  input.counter2TimesUs{trN} = NaN;
  input.counter2Values{trN} = NaN;
end



%% Block to acquire all sp/itiCounter times and tValues
try
  input.spCounter1TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter1', 'all', [], 1);
  input.spCounter1Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter1', 'all', 1) ;
catch
  input.spCounter1TimesUs{trN} = NaN;
  input.spCounter1Values{trN} = NaN;
end

try
  input.spCounter2TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter2', 'all', [], 1);
  input.spCounter2Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter2', 'all', 1) ;
catch
  input.spCounter2TimesUs{trN} = NaN;
  input.spCounter2Values{trN} = NaN;
end

try
  input.spCounter3TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter3', 'all', [], 1);
  input.spCounter3Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter3', 'all', 1) ;
catch
  input.spCounter3TimesUs{trN} = NaN;
  input.spCounter3Values{trN} = NaN;
end

try
  input.spCounter4TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter4', 'all', [], 1);
  input.spCounter4Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter4', 'all', 1) ;
catch
  input.spCounter4TimesUs{trN} = NaN;
  input.spCounter4Values{trN} = NaN;
end

try
  input.spCounter5TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter5', 'all', [], 1);
  input.spCounter5Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter5', 'all', 1) ;
catch
  input.spCounter5TimesUs{trN} = NaN;
  input.spCounter5Values{trN} = NaN;
end

try
  input.spCounter6TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter6', 'all', [], 1);
  input.spCounter6Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter6', 'all', 1) ;
catch
  input.spCounter6TimesUs{trN} = NaN;
  input.spCounter6Values{trN} = NaN;
end

try
  input.spCounter7TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter7', 'all', [], 1);
  input.spCounter7Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter7', 'all', 1) ;
catch
  input.spCounter7TimesUs{trN} = NaN;
  input.spCounter7Values{trN} = NaN;
end

try
  input.spCounter8TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter8', 'all', [], 1);
  input.spCounter8Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter8', 'all', 1) ;
catch
  input.spCounter8TimesUs{trN} = NaN;
  input.spCounter8Values{trN} = NaN;
end

try
  input.spCounter9TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter9', 'all', [], 1);
  input.spCounter9Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter9', 'all', 1) ;
catch
  input.spCounter9TimesUs{trN} = NaN;
  input.spCounter9Values{trN} = NaN;
end

try
  input.spCounter10TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'spCounter10', 'all', [], 1);
  input.spCounter10Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'spCounter10', 'all', 1) ;
catch
  input.spCounter10TimesUs{trN} = NaN;
  input.spCounter10Values{trN} = NaN;
end

try
  input.itiCounter1TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter1', 'all', [], 1);
  input.itiCounter1Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter1', 'all', 1) ;
catch
  input.itiCounter1TimesUs{trN} = NaN;
  input.itiCounter1Values{trN} = NaN;
end

try
  input.itiCounter2TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter2', 'all', [], 1);
  input.itiCounter2Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter2', 'all', 1) ;
catch
  input.itiCounter2TimesUs{trN} = NaN;
  input.itiCounter2Values{trN} = NaN;
end

try
  input.itiCounter3TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter3', 'all', [], 1);
  input.itiCounter3Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter3', 'all', 1) ;
catch
  input.itiCounter3TimesUs{trN} = NaN;
  input.itiCounter3Values{trN} = NaN;
end

try
  input.itiCounter4TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter4', 'all', [], 1);
  input.itiCounter4Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter4', 'all', 1) ;
catch
  input.itiCounter4TimesUs{trN} = NaN;
  input.itiCounter4Values{trN} = NaN;
end

try
  input.itiCounter5TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter5', 'all', [], 1);
  input.itiCounter5Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter5', 'all', 1) ;
catch
  input.itiCounter5TimesUs{trN} = NaN;
  input.itiCounter5Values{trN} = NaN;
end

try
  input.itiCounter6TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter6', 'all', [], 1);
  input.itiCounter6Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter6', 'all', 1) ;
catch
  input.itiCounter6TimesUs{trN} = NaN;
  input.itiCounter6Values{trN} = NaN;
end

try
  input.itiCounter7TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter7', 'all', [], 1);
  input.itiCounter7Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter7', 'all', 1) ;
catch
  input.itiCounter7TimesUs{trN} = NaN;
  input.itiCounter7Values{trN} = NaN;
end

try
  input.itiCounter8TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter8', 'all', [], 1);
  input.itiCounter8Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter8', 'all', 1) ;
catch
  input.itiCounter8TimesUs{trN} = NaN;
  input.itiCounter8Values{trN} = NaN;
end

try
  input.itiCounter9TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter9', 'all', [], 1);
  input.itiCounter9Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter9', 'all', 1) ;
catch
  input.itiCounter9TimesUs{trN} = NaN;
  input.itiCounter9Values{trN} = NaN;
end

try
  input.itiCounter10TimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'itiCounter10', 'all', [], 1);
  input.itiCounter10Values{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'itiCounter10', 'all', 1) ;
catch
  input.itiCounter10TimesUs{trN} = NaN;
  input.itiCounter10Values{trN} = NaN;
end


input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');
 trStr = sprintf('Trial %1.0f Completed', trN);
disp(trStr)

%% Run Plotting Function
input.holdStartsMs = input.tThisTrialStartTimeMs;
%input = exptRunSubfunctions(ds, input, { @plotVisStim });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Returns calculated and passed variables for the next trial.
retval = input;
return

