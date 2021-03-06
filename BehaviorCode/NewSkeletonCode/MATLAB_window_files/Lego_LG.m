function [retval] = Lego(data_struct, input)
% Main Matlab interface function for Lego-based experiments.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions
% Created 140612 by Andrew McKinney

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');

varsOneValEachTrial = { ...
    'tBlock2TrialNumber', ...
    'tLeftTrial', ...
    'tStimProbAvgLeft', ...
    'tLeftBias', ...
    'tLastTrialWasLeft', ...
    'tGratingEccentricityStart', ...
    'tGratingContrast', ...
    'dGratingContrast', ...
    'dGratingContrastDiff', ...
    'tGratingDiameterDeg', ...
    'dGratingDiameterDeg', ...
    'dGratingDiameterDiff', ...
    'tGratingStartingPhaseDeg', ...
    'tDecisionTimeMs', ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tConsecCorrects', ...
    'tConsecErrors', ...
    'tNTrialsCompleted', ...
    'tNBlockLeftTrsCompleted', ...
    'tNBlockRightTrsCompleted', ...
    'tNStimAccepted', ...
    'tRewardTimeUs', ...
    'gratingDirectionDeg', ...
    'tLaserPowerMw', ...
    'tLaserDoLinearRamp', ...
    'tLaserDoPulseTrain', ...
    'tLaserRampLengthMs', ...
    'tLaserPeakMaxMw', ...
    'tLaserBaselinePowerMw', ...
    'tLaserPeakStepsPerOctave', ...
    'tTrialLaserPowerMw', ...
    'tTrialLaserOnTimeMs', ...
    'tTrialLaserOffTimeMs', ...
    'robotGoLeft', ...
    'stimTimestampMs', ...
    'qStimOn', ...
    'qTrialStart', ...
    'qStartReact', ...
    'isNoGo', ...
    'didNoGo', ...
    'inPreBlock', ...
    'tPuffDurationMs', ...
    'tLeftResponse', ...
    'tRightResponse', ...
    'tStationaryPeriodMs', ...
    'tItiTimeMs', ...
    'tGratingDirectionStart', ...
    'tFeedbackMotionSensitivity', ...
    'adapterTimestampMs', ...
    'aGratingContrast', ...
    'aGratingDirectionDeg', ...
    'aGratingPhaseDeg', ...
    'tAGratingDirectionDeg', ...
    'tDoFeedbackMotion', ...
    'tDoMask',...
    'tMaskDirectionDeg', ...
    'tMaskDirectionStart', ...
    'tMaskStartingPhaseDeg', ...
    'dMaskStartingPhaseDeg', ...
    'tMaskContrast', ...
    'dMaskContrast', ...
    'tPlaidDirectionDeg', ...
    'tPlaidDirectionStart', ...
    'tGratingSpatialFreqCPD', ...

  };

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% process reaction times for this trial
stimOnUs = mwGetEventTime(eventsTrial, ds.event_codec, 'stimulusOn', 1);
%try
%    input.quadStampsUs{trN} = [0 ((mwGetEventTime(eventsTrial, ds.event_codec, 'quadrature', 'all', [], 1))-stimOnUs)/1000];
%    input.quadValues{trN} = [0 (mwGetEventValue(eventsTrial, ds.event_codec, 'quadrature', 'all', 1)) - input.tQuadrature{trN}];
%catch
%    input.quadStampsUs{trN} = NaN;
%    input.quadValues{trN} = NaN;
%end

input.allAdaptGratingDirectionDeg{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'aGratingDirectionDeg', 'all', 1);
input.allAdaptGratingPhaseDeg{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'aGratingPhaseDeg', 'all', 1);

if input.tLeftTrial{trN}==1,
    input.leftGratingContrast{trN} = input.tGratingContrast{trN};
    input.rightGratingContrast{trN} = input.dGratingContrast{trN};
    input.leftGratingDiameterDeg{trN} = input.tGratingDiameterDeg{trN};
    input.rightGratingDiameterDeg{trN} = input.dGratingDiameterDeg{trN};
else
    input.leftGratingContrast{trN} = input.dGratingContrast{trN};
    input.rightGratingContrast{trN} = input.tGratingContrast{trN};
    input.leftGratingDiameterDeg{trN} = input.dGratingDiameterDeg{trN};
    input.rightGratingDiameterDeg{trN} = input.tGratingDiameterDeg{trN};
end

try
    input.quadratureTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'quadrature', 'all', [], 1);
    input.quadratureValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'quadrature', 'all', 1) ;
catch
    input.quadratureTimesUs{trN} = NaN;
    input.quadratureValues{trN} = NaN;
end

try
    input.tGratingDirectionDeg{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'tGratingDirectionDeg', 'all', 1) ;
    input.tGratingEccentricityDeg{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'tGratingEccentricityDeg', 'all', 1) ;
    input.dGratingEccentricityDeg{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'dGratingEccentricityDeg', 'all', 1) ;
catch
    input.tGratingDirectionDeg{trN} = NaN;
    input.tGratingEccentricityDeg{trN} = NaN;
    input.dGratingEccentricityDeg{trN} = NaN;
end

%% disp status
if input.tLeftTrial{trN}==1,
    tLeftStr = 'Left Trial';
else
    tLeftStr = 'Right Trial';
end
outcomeStr = input.trialOutcomeCell{trN};
outcomeStr = strcat(upper(outcomeStr(1)), outcomeStr(2:end));
decisionTime = input.tDecisionTimeMs{trN};
if isempty(input.tRewardTimeUs{trN}),
    rewS = 0;
else
    rewS = input.tRewardTimeUs{trN}/1000;
end

if input.doContrastDiscrim
fprintf(1,'Con: T=%0.2f, D=%0.2f, %s, %s, RT: %0.0f ms, B%1.0f, Rew: %2.0f\n ', input.tGratingContrast{trN}, input.dGratingContrast{trN}, tLeftStr, outcomeStr, decisionTime, input.tBlock2TrialNumber{trN}+1, rewS)
elseif input.doSizeDiscrim
fprintf(1,'Size: T=%2.0f, D=%2.0f, %s, %s, RT: %0.0f ms, B%1.0f, Rew: %2.0f\n ', input.tGratingDiameterDeg{trN}, input.dGratingDiameterDeg{trN}, tLeftStr, outcomeStr, decisionTime, input.tBlock2TrialNumber{trN}+1, rewS)
elseif input.doOriDiscrim
fprintf(1,'Ori: T=%2.0f, %s, %s, RT: %0.0f ms, B%1.0f, Rew: %2.0f\n ', abs(input.gratingTargetDirection - input.tGratingDirectionStart{trN}), tLeftStr, outcomeStr, decisionTime, input.tBlock2TrialNumber{trN}+1, rewS)
end
% input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
%                               '~/Documents/MWorks/Data', ...
%                               mat2str(input.subjectNum), input.saveTime);

%itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
%fprintf(1,'Hold %d, req %d, react %d, %s%s %s- %d rew %dms\n', ...
%        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
%        itiStr, ...
%        stimStr, block2Str, ...
%        nJ, round(juiceD));

%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotLego });

save(input.savedDataName, 'input');

%% save variables for next trial
retval = input;
%input = exptSaveMatlabState(data_struct, input);

return
