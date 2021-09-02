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
    'tItiWaitTimeMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tTempStim1', ...
    'tTempStim2', ...
    'tTempMask1', ...
    'tTempMask2', ...
    'tStimulusNumber', ...
    'tStimulusNumber2', ...
    'tMaskNumber', ...
    'tMaskNumber2', ...
    'tStimOneGratingOnTimeMs', ...
    'tStimTwoGratingOnTimeMs', ...
    'tDoMask',...
    'tISITimeMs', ...
    'tStimOneDoVisualStim', ...
    'tStimOneDoAuditoryStim', ...
    'tStimOneGratingContrast', ...
    'tStimOneGratingDirectionDeg', ...
    'tStimOneGratingElevationDeg', ...
    'tStimOneGratingAzimuthDeg', ...
    'tStimOneGratingDiameterDeg', ...
    'tStimOneGratingSpatialFreqCPD', ...
    'tStimOneGratingTemporalFreqCPS', ...
    'tStimOneGratingPhaseDeg', ...
    'tStimOneSoundAmplitude', ...
    'tMaskOneGratingContrast', ...
    'tMaskOneGratingDirectionDeg', ...
    'tMaskOneGratingPhaseDeg', ...
    'tMaskOneGratingTemporalFreqCPS', ...
    'tStimTwoDoVisualStim', ...
    'tStimTwoDoAuditoryStim', ...
    'tStimTwoGratingContrast', ...
    'tStimTwoGratingDirectionDeg', ...
    'tStimTwoGratingElevationDeg', ...
    'tStimTwoGratingAzimuthDeg', ...
    'tStimTwoGratingDiameterDeg', ...
    'tStimTwoGratingSpatialFreqCPD', ...
    'tStimTwoGratingTemporalFreqCPS', ...
    'tStimTwoGratingPhaseDeg', ...
    'tStimTwoSoundAmplitude', ...
    'tMaskTwoGratingContrast', ...
    'tMaskTwoGratingDirectionDeg', ...
    'tMaskTwoGratingPhaseDeg', ...
    'tMaskTwoGratingTemporalFreqCPS', ...
    'tDoRotate', ...
    'tDoMaskRotate', ...
    'tDoMatrix', ...
    'tDoRandISITime', ...
    'tDoRandStimOnTime', ...
    'tDoRandCon', ...
    'tDoRandMaskCon', ...
    'tDoRandDir', ...
    'tDoRandMaskDir', ...
    'tDoRandPos', ...
    'tDoRandSF', ...
    'tDoRandTF', ...
    'tDoRandPhase', ...
    'tDoRandMaskPhase', ...
    'tDoRandMaskTF', ...
    'tDoRandDiameter',...
    'tDoRandSoundAmp',...
    'tTrialLaserPowerMw', ...
    'tTrialLaserOnTimeMs', ...
    'tTrialLaserOffTimeMs', ...
    'tBlock2TrialNumber', ...
    'tNStimAccepted', ...
    'tItiWaitFrames', ...
    'nStimOneFramesOn', ...
    'nStimTwoFramesOn', ...
    'nFramesISI', ...
    'cItiStart',...
    'cStimOneOn',...
    'cStimTwoOn',...
    'cStimOneOff',...
    'cStimTwoOff',...
    'mwStimOneOnMs',...
    'mwStimTwoOnMs',...
    'mwStimOneOffMs',...
    'mwStimTwoOffMs',...
    'tDoRandSoundFreq',...
    'tStimOneSound1Amplitude',...
    'tStimOneSound2Amplitude',...
    'tStimOneSound3Amplitude',...
    'tStimOneSound4Amplitude',...
    'tStimOneSound5Amplitude',...
    'tStimTwoSound1Amplitude',...
    'tStimTwoSound2Amplitude',...
    'tStimTwoSound3Amplitude',...
    'tStimTwoSound4Amplitude',...
    'tStimTwoSound5Amplitude',...
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
    'doMatrix', ...
    'doMask', ...
    'doTwoStimTogether', ...
    'fractMaskTrials', ...
    'doRandISITime', ...
    'doRandStimOnTime', ...
    'doRandCon', ...
    'doZeroCon',...
    'zeroConPer80',...
    'doRandDir', ...
    'doRandMaskCon', ...
    'doRandMaskDir', ...
    'doRandMaskPhase', ...
    'doRandMaskTF', ...
    'doRandPos', ...
    'doRandDiameter', ...
    'doRandSF', ...
    'doRandTF', ...
    'doRandPhase', ...
    'doRotate',...
    'doMaskRotate',...
    'doRandSoundAmp', ...
    'doWheelSpeed', ...
    'doBlock2', ...
    'doBlock2SeparateOdds', ...
    'block2TrPer80', ...
    'block2MatchB1VisStim', ...
    'block2MatchB1AudStim', ...
    'block2ISITimeMs', ...
    'block2DoMatrix', ...
    'block2DoRandISITime', ...
    'block2DoRandStimOnTime', ...
    'block2DoRandCon', ...
    'block2DoRandDir', ...
    'block2DoRandPos', ...
    'block2DoRandDiameter', ...
    'block2DoRandSF', ...
    'block2DoRandTF', ...
    'block2DoRandPhase', ...
    'block2DoRotate', ...
    'block2DoMask', ...
    'block2DoMaskRotate', ...
    'block2DoRandSoundAmp', ...
    'frameRateHz', ...
    'gratingRotateSpeedDPS',...
    'itiTimeMs', ...
    'doExtendItiOnShortPrevTrial', ...
    'stimOneDoVisualStim', ...
    'stimOneDoAuditoryStim', ...
    'stimOneGratingOnTimeMs', ...
    'stimOneGratingContrast', ...
    'stimOneGratingDirectionDeg', ...
    'stimOneGratingElevationDeg', ...
    'stimOneGratingAzimuthDeg', ...
    'stimOneGratingDiameterDeg', ...
    'stimOneGratingSpatialFreqCPD', ...
    'stimOneGratingTemporalFreqCPS', ...
    'stimOneGratingPhaseDeg', ...
    'stimOneSoundAmplitude', ...
    'maskOneGratingContrast', ...
    'maskOneGratingDirectionDeg', ...
    'maskOneGratingPhaseDeg', ...
    'maskOneGratingTemporalFreqCPS', ...
    'stimTwoDoVisualStim', ...
    'stimTwoDoAuditoryStim', ...
    'stimTwoGratingOnTimeMs', ...
    'stimTwoGratingContrast', ...
    'stimTwoGratingDirectionDeg', ...
    'stimTwoGratingElevationDeg', ...
    'stimTwoGratingAzimuthDeg', ...
    'stimTwoGratingDiameterDeg', ...
    'stimTwoGratingSpatialFreqCPD', ...
    'stimTwoGratingTemporalFreqCPS', ...
    'stimTwoGratingPhaseDeg', ...
    'stimTwoSoundAmplitude', ...
    'maskTwoGratingContrast', ...
    'maskTwoGratingDirectionDeg', ...
    'maskTwoGratingPhaseDeg', ...
    'maskTwoGratingTemporalFreqCPS', ...
    'matchStimOneParameters', ...
    'block2DoTrialLaser', ...
    'block2TrialLaserPowerMw', ...
    'block2TrialLaserOnTimeMs', ...
    'block2TrialLaserOffTimeMs', ...
    'block2StimOneDoVisualStim', ...
    'block2StimOneDoAuditoryStim', ...
    'block2StimOneGratingOnTimeMs', ...
    'block2StimOneGratingContrast', ...
    'block2StimOneGratingDirectionDeg', ...
    'block2StimOneGratingElevationDeg', ...
    'block2StimOneGratingAzimuthDeg', ...
    'block2StimOneGratingDiameterDeg', ...
    'block2StimOneGratingSpatialFreqCPD', ...
    'block2StimOneGratingTemporalFreqCPS', ...
    'block2StimOneGratingPhaseDeg', ...
    'block2StimOneSoundAmplitude', ...
    'block2MaskOneGratingContrast', ...
    'block2MaskOneGratingDirectionDeg', ...
    'block2StimTwoDoVisualStim', ...
    'block2StimTwoDoAuditoryStim', ...
    'block2StimTwoGratingOnTimeMs', ...
    'block2StimTwoGratingContrast', ...
    'block2StimTwoGratingDirectionDeg', ...
    'block2StimTwoGratingElevationDeg', ...
    'block2StimTwoGratingAzimuthDeg', ...
    'block2StimTwoGratingDiameterDeg', ...
    'block2StimTwoGratingSpatialFreqCPD', ...
    'block2StimTwoGratingTemporalFreqCPS', ...
    'block2StimTwoGratingPhaseDeg', ...
    'block2StimTwoSoundAmplitude', ...
    'block2MaskTwoGratingContrast', ...
    'block2MaskTwoGratingDirectionDeg', ...
    'block2MatchStimOneParameters',...
    'doRandSoundFreq',...
    'block2DoRandSoundFreq',...
    'stimOneSound1Amplitude',...
    'stimOneSound2Amplitude',...
    'stimOneSound3Amplitude',...
    'stimOneSound4Amplitude',...
    'stimOneSound5Amplitude',...
    'stimTwoSound1Amplitude',...
    'stimTwoSound2Amplitude',...
    'stimTwoSound3Amplitude',...
    'stimTwoSound4Amplitude',...
    'stimTwoSound5Amplitude',...
    'block2StimOneSound1Amplitude',...
    'block2StimOneSound2Amplitude',...
    'block2StimOneSound3Amplitude',...
    'block2StimOneSound4Amplitude',...
    'block2StimOneSound5Amplitude',...
    'block2StimTwoSound1Amplitude',...
    'block2StimTwoSound2Amplitude',...
    'block2StimTwoSound3Amplitude',...
    'block2StimTwoSound4Amplitude',...
    'block2StimTwoSound5Amplitude',...
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

try
    input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
    input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
catch
    input.wheelSpeedTimesUs{trN} = NaN;
    input.wheelSpeedValues{trN} = NaN;
end

%% disp status

trStr = ['Trial ' num2str(trN) '- '];
if input.tStimOneGratingContrast{trN} > 0
    stim1Str = sprintf('Stim1: con %g%%, dir %g%, %g% ms', input.tStimOneGratingContrast{trN}.*100, input.tStimOneGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    stim1Str = sprintf('Stim1: none');
end
if input.tMaskOneGratingContrast{trN} > 0
    mask1Str = sprintf('Mask1: con %g%%, dir %g%, %g% ms', input.tMaskOneGratingContrast{trN}*100, input.tMaskOneGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    mask1Str = sprintf('Mask1: none');
end
if input.tStimTwoGratingContrast{trN} > 0
    stim2Str = sprintf('Stim2: con %g%%, dir %g%, %g% ms', input.tStimTwoGratingContrast{trN}.*100, input.tStimTwoGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    stim2Str = sprintf('Stim2: none');
end
if input.tMaskTwoGratingContrast{trN} > 0
    mask2Str = sprintf('Mask2: con %g%%, dir %g%, %g% ms', input.tMaskTwoGratingContrast{trN}*100, input.tMaskTwoGratingDirectionDeg{trN}, input.tStimTwoGratingOnTimeMs{trN});
else
    mask2Str = sprintf('Mask2: none');
end
isiStr = sprintf('ISI: %g%ms', input.tISITimeMs{trN});

if ~input.doBlock2
  block2Str = '';
else
  block2Str = sprintf(' b2tr %d ', input.tBlock2TrialNumber{trN});
end
itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
fprintf(1,'%s %s; %s; %s; %s %s\n', ...
    trStr, stim1Str, stim2Str, mask1Str, mask2Str, isiStr, block2Str);

%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotTwoStim_2P_Frames });


input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');

%% save variables for next trial
retval = input;

return


