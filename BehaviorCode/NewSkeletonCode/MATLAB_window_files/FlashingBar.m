function [retval] = MovingBar(data_struct, input)
%Main MATLAB interface function for VisStimRet.xml.
%
if nargin < 2,
input = [];
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');

% Vars set on each trial
oneValEachTrialNames = { ...
    'counter', ...
    'tTrialsDoneSinceStart' , ...
    'tXStartPos',...
    'tYStartPos',...
    'tXSizeDeg',...
    'tYSizeDeg',...
    'tXDirection',...
    'tYDirection',...
    'tBarOrientationDeg',...
    'tStimulusNumber',...
    'tNStimAccepted',...
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
input.constList = { 'subjectNum',...
    'backgroundLuminance',...
    'barWidth',...
    'monitorXDeg',...
    'monitorYDeg',...
    'barOrientationDeg',...
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

%% Counter/Frames Synchronization Section
try
input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
input.counterTimesUs{trN} = NaN;
input.counterValues{trN} = NaN;
end

%% Running events Section
try
input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
catch
input.wheelSpeedTimesUs{trN} = NaN;
input.wheelSpeedValues{trN} = NaN;
end


input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');


if input.tBarOrientationDeg{trN} == 0
dir_str = 'right';
pos_idx = num2str(input.tXStartPos{trN});
elseif input.tBarOrientationDeg{trN} == 90
dir_str = 'up';
pos_idx = num2str(input.tYStartPos{trN});
end
disp(sprintf('Trial %d: Ori %d, Direction %s,Position %s ', ...
             trN, input.tBarOrientationDeg{trN}, dir_str, pos_idx));

%% Run Plotting Function
input = exptRunSubfunctions(ds, input, { @plotFlashingBar });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%trN = input.trialSinceReset;

%stimStr = [stimStr sprintf('Contrast %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];

%% Returns calculated and passed variables for the next trial.
retval = input;
return

