function [retval] = RewardDelay(data_struct, input)
% Online Matlab processing for RewardDelay.xml
%
% Created June 25, 2013 by Andrew McKinney
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2, input = []; end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');

varsOneValEachTrial = { ...
    ...%'holdStartsMs', ...
    'tTotalReqHoldTimeMs', ...
    'tTotalRewardTimeUs', ...
    'tRandReqHoldTimeMs', ...
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

%% process reaction times for this trial
stimOnUs = mwGetEventTime(eventsTrial, ds.event_codec, 'stimulusOn', 1);

%trialStartUs = mwGetEventTime(eventsTrial, ds.event_codec, 'trialStart', 1);
leverDownUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 1, 1);
leverUpUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 2, 0);

%reactTimeMs = (leverUpUs - stimOnUs) / 1000;
totalHoldTimeMs = input.tTotalReqHoldTimeMs{trN};
holdTimeMs = (leverUpUs - leverDownUs) / 1000;
reactTimeMs = (holdTimeMs - totalHoldTimeMs);

% add to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;


%% disp status
juiceAmtsMs = input.juiceTimesMsCell{trN};
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);


if input.itiExtraRandTimeMs > 0
  itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
else
  itiStr = '';
end
fprintf(1,'Hold %d, req %d, react %d, %s%s %s- %d rew %dms\n', ...
        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
        itiStr, ...
        nJ, round(juiceD));

%% Run Reward Delay Plotting Function
input = exptRunSubfunctions(ds, input, { @plotRewardDelay });



%% Save Variables for next trial
retval = input;

return


