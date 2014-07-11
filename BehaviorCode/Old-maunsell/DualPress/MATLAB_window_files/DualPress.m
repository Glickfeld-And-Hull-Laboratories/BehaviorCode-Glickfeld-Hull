function [retval] = DualPress(data_struct, input)
% Main matlab online function for DualPress
%
% MH 130107: refactored into new expt* functions.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2, input = []; end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/Old-maunsell/MatlabSharedCode');

dualPressVarsOneValEachTrial = { ...
    'tGratingContrast', ... 
    'tMarkovOdds', ... 
    'tGratingDirectionDeg', ... 
    'tGratingElevationDeg', ... 
    'tGratingAzimuthDeg', ...
    'tGratingHeightDeg', ...
    'tGratingWidthDeg', ...
    'tGratingSpatialFreqCPD', ...
    'tGratingSpeedDPS', ...
    'tGratingDurationMs', ...
    'tLeftTrial', ...
    'tFirstReactReleaseIsLeft', ...
    'tThisTrialStartTimeMs', ...
    'tStartTrialWaitForPressTimeMs', ...
    'tRewardRunningLeftBias', ...
    'tTotalRewardTimeUs', ...
    'tTotalReqHoldTimeMs', ...
    'FIO1', ...
    'FIO4', ...
    'pressTimestampMs', ...

    };

... %'leverLatencyRPosMs', ...
... %'responseSideIx', ...

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, dualPressVarsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;


%% process reaction times for this trial

% backward compat
totalHoldTimeMs = input.tTotalReqHoldTimeMs{trN};

trialStartUs = mwGetEventTime(eventsTrial, ds.event_codec, 'trialStart', 'all', 1, false);
leverDownUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 'all', 1, false);
leverUpUs = mwGetEventTime(eventsTrial, ds.event_codec, 'leverResult', 2, 0, false);  % occ 1 happens always at trial start
assert(length(trialStartUs) == 1);
assert(length(leverDownUs) == 1);

holdTimeMs = (leverUpUs - leverDownUs) / 1000;
reactTimeMs = (holdTimeMs - totalHoldTimeMs);

% add to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;
input.reqHoldTimeMs{trN} = totalHoldTimeMs;

juiceAmtsMs = input.juiceTimesMsCell{trN};


% contrast
input.gratingContrast{trN} = input.tGratingContrast{trN};

% laser power does not work now
% $$$ if input.doRightLaserStim == 1
% $$$     % laser power values  
% $$$     tLaserPowerMw = mwGetEventValue(eventsTrial, ds.event_codec, ...
% $$$                                     'tTrialLaserPowerMw', 'all');
% $$$     input.laserPowerMw{thisTrialN} = tLaserPowerMw(end);
% $$$     %disp(sprintf('Power was %3.2fmW', input.laserPowerMw{thisTrialN}));
% $$$ end

%% parse response side
input.lever1State{trN}=mwGetEventValue(eventsTrial, ds.event_codec, 'FIO1', 'all', 'ignoreMissing');
input.lever2State{trN}=mwGetEventValue(eventsTrial, ds.event_codec, 'FIO4', 'all', 'ignoreMissing');
input.lever1Ts{trN}=mwGetEventTime(eventsTrial, ds.event_codec, 'FIO1', 'all', [], 'ignoreMissing');
input.lever2Ts{trN}=mwGetEventTime(eventsTrial, ds.event_codec, 'FIO4', 'all', [], 'ignoreMissing');

holdStartMs = input.pressTimestampMs{trN};
lever1Ts = input.lever1Ts{trN};
lever2Ts = input.lever2Ts{trN};
tLeftTrial = input.tLeftTrial{trN};

firstPossReleaseMs = double(holdStartMs + input.reqHoldTimeMs{trN});
if strcmp(input.trialOutcomeCell{trN}, 'early')
    firstPossReleaseMs = double(holdStartMs + holdTimeMs - 5);  % count from half-a-poll before early event
end

first1N = find(lever1Ts > 1000*firstPossReleaseMs, 1, 'first');
first2N = find(lever2Ts > 1000*firstPossReleaseMs, 1, 'first');
leverRReactMs = double(lever1Ts(first1N))/1000 - firstPossReleaseMs;
leverLReactMs = double(lever2Ts(first2N))/1000 - firstPossReleaseMs;


% deal with case where only one lever was actually released
if isempty(first1N)
    leverRReactMs = NaN;
end
if isempty(first2N)
    leverLReactMs = NaN;
end

% find latency and save
leverLatencyRPosMs = leverLReactMs - leverRReactMs;
input.leverLatencyRPosMs{trN} = leverLatencyRPosMs;

%% display text string
if strcmp(input.trialOutcomeCell{trN}, 'early')
    outcomeStr = 'E';
elseif strcmp(input.trialOutcomeCell{trN}, 'ignore')
    outcomeStr = 'M';
elseif strcmp(input.trialOutcomeCell{trN}, 'dualrelease') 
    outcomeStr = 'D';
elseif strcmp(input.trialOutcomeCell{trN}, 'success') && tLeftTrial == 1 
    outcomeStr = 'L';
elseif strcmp(input.trialOutcomeCell{trN}, 'success') && tLeftTrial == 0
    outcomeStr =  'R';
elseif strcmp(input.trialOutcomeCell{trN}, 'incorrect') && tLeftTrial == 1
    outcomeStr = 'R';
elseif strcmp(input.trialOutcomeCell{trN}, 'incorrect') && tLeftTrial == 0
    outcomeStr = 'L';
end

% disp status
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
stimStr = '';
stimStr = [stimStr sprintf('ctst %g%%', chop(input.gratingContrast{trN},2)*100)];
% if input.doRightLaserStim
%     stimStr = [stimStr sprintf('pow %gmW', ...
%                                chop(input.laserPowerMw{trN}, 2))];
% end
block2Str = '-';
if tLeftTrial == 1
    stimSideStr = 'L';
elseif tLeftTrial == 0
    stimSideStr = 'R';
else
    error('tLeftTrial should be defined every trial');
end

fprintf(1,'H %d - req %d = %dms, %s St resp %s %s off+r %dms %d rew %dms\n', ...
        round(holdTimeMs), round(totalHoldTimeMs), round(reactTimeMs), ...
        stimStr, stimSideStr, outcomeStr, round(leverLatencyRPosMs), ...
        nJ, round(juiceD));



%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotDualPressHist });



%% save variables for next trial
retval = input;


return

