function [retval] = RewardDelay(data_struct, input)
% Online Matlab processing for RewardDelay.xml
% Code for Hull Lab - Duke Neurobiology
% Created June 25, 2013 by Andrew McKinney
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2, input = []; end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');

varsOneValEachTrial = { ...
    ...%'holdStartsMs', ...
    'tTotalReqHoldTimeMs', ...
    'tTotalRewardTimeUs', ...
    'tFakeMouseReactMs', ...
    'tFakeMousePressMs', ...
    'tStartTrialWaitForPressTimeMs', ...
    'tConsecErrors', ...
    'tConsecTimeoutStartTime', ...
    'tNRewards', ...
    'tInterRewardIntervalMs', ...
    'tRewardAddPerMsHoldUs', ...
    'tStimTurnedOn', ...
    'tItiWaitTimeMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs'
};

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



trN = input.trialSinceReset;

%% Process reaction and hold times for this trial
stimOnUs = mwGetEventTime(eventsTrial, ds.event_codec, 'stimulusOn', 1);
leverDownUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 1, 1);
leverUpUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 2, 0);

totalHoldTimeMs = input.tTotalReqHoldTimeMs{trN};
holdTimeMs = (leverUpUs - leverDownUs) / 1000;
reactTimeMs = (holdTimeMs - totalHoldTimeMs);

% Add calculated values to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;

%% MATLAB Interface Window display updating
juiceAmtsMs = input.juiceTimesMsCell{trN};
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);

fprintf(1,'ActualHold %d, Reqd %d, React %d, ITI %d - %d Reward Size %dms\n', ...
        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
        round(input.tItiWaitTimeMs{trN}), ...
        nJ, round(juiceD));

%% Run Reward Delay Plotting Function
input = exptRunSubfunctions(ds, input, { @plotRewardDelay });

%% Save Variables for next trial
retval = input;

return


