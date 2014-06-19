function [retval] = Lego(data_struct, input)
% Main Matlab interface function for Lego-based experiments.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions
% Created 140612 by Andrew McKinney

if nargin < 2,
    input = struct;
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');

varsOneValEachTrial = { ...
    'tLeftTrial', ...
    'tLastTrialWasLeft', ...
    'tGratingEccentricityDeg', ...
    'dGratingEccentricityDeg', ...
    'tGratingContrast', ...
    'dGratingContrast', ...
    'tQuadrature', ...
    'tRightReversal', ...
    'tLeftReversal', ...
    'tDecisionTimeMs', ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tConsecCorrects', ...
    'tConsecErrors', ...
    'tNTrialsCompleted', ...
  };

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% process reaction times for this trial
stimOnUs = mwGetEventTime(eventsTrial, ds.event_codec, 'stimulusOn', 1);
input.quadStampsUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'quadrature', 'all', [], 1);

%% disp status
if input.tLeftTrial{trN}==1,
    tLeftStr = 'Left Trial';
else
    tLeftStr = 'Right Trial';
end
outcomeStr = input.trialOutcomeCell{trN};
outcomeStr = strcat(upper(outcomeStr(1)), outcomeStr(2:end));

fprintf(1,'Target Contrast: %0.2f, Distractor Contrast: %0.2f, %s, %s\n ', input.tGratingContrast{trN}, input.dGratingContrast{trN}, tLeftStr, outcomeStr)


%itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
%fprintf(1,'Hold %d, req %d, react %d, %s%s %s- %d rew %dms\n', ...
%        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
%        itiStr, ...
%        stimStr, block2Str, ...
%        nJ, round(juiceD));

%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotLego });

%% save variables for next trial
retval = input;
input = exptSaveMatlabState(data_struct, input);

return
