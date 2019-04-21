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
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tTempStim', ...
    'tTempMask', ...
    'tStimulusNumber', ...
    'tMaskNumber', ...
    'tDoContingentAdapt', ...
    'tDoAsynchronousAdapt', ...
    'tTestGratingContrast', ...
    'tTestGratingDirectionDeg', ...
    'tMaskGratingContrast', ...
    'tMaskGratingDirectionDeg', ...
    'tAdaptGratingContrast', ...
    'tAdaptMaskGratingContrast', ...
    'tAdaptGratingDirectionDeg', ...
    'tAdaptMaskGratingDirectionDeg', ...
    'tAdaptTimeMs',...
    'tDoRandCon', ...
    'tDoRandDir', ...
    'tNStimAccepted', ...
    'tNMaskAccepted', ...
    'tItiTimeMs', ...
    'mwTestOnMs',...
    'mwTestOffMs',...
    'quadrature',...
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
    'speedIntervalMS', ...
    'doRandCon', ...
    'doRandDir', ...
    'doWheelSpeed', ...
    'itiTimeMs', ...
    'initialAdaptTimeMs',...
    'topUpAdaptTimeMs',...
    'adaptStimOnTimeMs', ...
    'adaptGratingContrast',...
    'adaptGratingDirectionDeg',...
    'adaptMaskGratingDirectionDegDiff',...
    'trialsPerAdaptBlock',...
    'testStimOnTimeMs', ...
    'testGratingContrast', ...
    'testGratingStepLog', ...
    'testGratingStepN', ...
    'maskGratingContrast', ...
    'maskGratingStepLog', ...
    'maskGratingStepN', ...
    'testGratingDirectionDeg', ...
    'testGratingDirectionDegStep', ...
    'testGratingDirectionDegStepN', ...
    'maskGratingDirectionDeg', ...
    'maskGratingDirectionDegStep', ...
    'maskGratingDirectionDegStepN', ...
    'gratingElevationDeg', ...
    'gratingAzimuthDeg', ...
    'gratingDiameterDeg', ...
    'gratingSpatialFreqCPD', ...
    'gratingTemporalFreqCPS', ...
    'gratingPhaseDeg', ...
    'sendElectrophysSync',...
    'doPhotodiode'};

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

try
    input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
    input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
catch
    input.wheelSpeedTimesUs{trN} = NaN;
    input.wheelSpeedValues{trN} = NaN;
end

%% disp status
if input.tStimOneGratingContrast{trN} > 0
    stim1Str = sprintf('Stim1: con %g%%, dir %g%, %g% ms', input.tStimOneGratingContrast{trN}.*100, input.tStimOneGratingDirectionDeg{trN}, input.tStimOneGratingOnTimeMs{trN});
else
    stim1Str = sprintf('');
end
if input.tMaskOneGratingContrast{trN} > 0
    mask1Str = sprintf('Mask1: con %g%%, dir %g%, %g% ms', input.tMaskOneGratingContrast{trN}*100, input.tMaskOneGratingDirectionDeg{trN}, input.tStimOneGratingOnTimeMs{trN});
else
    mask1Str = sprintf('');
end
if input.tStimTwoGratingContrast{trN} > 0
    stim2Str = sprintf('Stim1: con %g%%, dir %g%, %g% ms', input.tStimTwoGratingContrast{trN}.*100, input.tStimTwoGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    stim2Str = sprintf('');
end
if input.tMaskTwoGratingContrast{trN} > 0
    mask2Str = sprintf('Mask1: con %g%%, dir %g%, %g% ms', input.tMaskTwoGratingContrast{trN}*100, input.tMaskTwoGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    mask2Str = sprintf('');
end
isiStr = sprintf('ISI: %g%ms', input.tISITimeMs{trN});

if ~input.doBlock2
  block2Str = '';
else
  block2Str = sprintf(' b2tr %d ', input.tBlock2TrialNumber{trN});
end
itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
fprintf(1,'%s %s %s %s %s %s\n', ...
    stim1Str, stim2Str, mask1Str, mask2Str, isiStr, block2Str);


%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotTwoStim });



%% save variables for next trial
retval = input;

return


