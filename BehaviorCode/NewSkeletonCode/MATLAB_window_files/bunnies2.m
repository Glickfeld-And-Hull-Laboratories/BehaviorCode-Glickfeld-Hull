function [retval] = HoldAndDetectConstant(data_struct, input)
% Main matlab online function for HADC8
%
%  MH 100115: created
%  MH 130107: refactored into new expt* functions.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2,
    input = [];
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');

oneValEachTrialNames = {...
    'ntrials', ...
    'cItiStart', ...
    'cStimOneOn', ...
    'cStimOneOff', ...
    'cStimTwoOn', ...
    'cStimTwoOff', ...
    'tisiTimeMs', ...
    'tstimOne', ...
    'tstimTwo', ...
    'tisiTimeFrames', ...
    'titiTimeFrames', ...
    'tstimOnTimeFrames', ...
};

events = ds.events;
eventsTrial = events;

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
input.constList = {'subjectNum',...
    'frameRateHz',...
    'index', ...
    'xSize', ...
    'ySize', ...
    'xPosition', ...
    'yPosition', ...
    'isiTimeMs', ...
    'randIsiTimeStart', ...
    'randIsiTimeInterval', ...
    'itiTimeMs', ...
    'stimOnTimeMs', ...
    'stimOne', ...
    'stimTwo', ...
    'doRandIsiTime', ...
    'doRepeatStimOne', ...
    'doRepeatStimTwo', ...
    'doRandStimOne', ...
    'doRandStimTwo', ...
    'doSameStims', ...
    'DoWheelSpeed', ...
};

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% Counter/Frames Synchronization Section
try
    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

% try
    % input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
   % input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
% catch
   % input.wheelSpeedTimesUs{trN} = NaN;
   % input.wheelSpeedValues{trN} = NaN;
% end

%% disp status

trStr = ['Trial ' num2str(trN) '- '];


%% run subfunctions
%input.changedStr(trN)
input = exptRunSubfunctions(ds, input, { @plotBunnies2 });

input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');

%% save variables for next trial
retval = input;

return


