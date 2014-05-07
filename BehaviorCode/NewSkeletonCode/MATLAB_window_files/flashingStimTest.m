function [retval] = HoldAndDetectConstant(data_struct, input)
% Main matlab online function for HADC8
%
%  MH 100115: created
%  MH 130107: refactored into new expt* functions.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2, input = []; end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/Old-maunsell/MatlabSharedCode');

varsOneValEachTrial = {...
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
    'tLastTrialStartTimeMs', ...
    'tTempStimOdds', ...
    'tStimulusNumber', ...
    'tSvStimNumber', ...
    'tLaserPowerMw', ...
    'tLaserDoLinearRamp', ...
    'tLaserDoPulseTrain', ...
    'tLaserRampLengthMs', ...
    'tLaserPeakMaxMw', ...
    'tLaserBaselinePowerMw', ...
    'tGratingMaxContrastStep', ...
    'tGratingContrastStepsPerOctave', ...
    'tGratingMaxDirectionStepDeg', ...
    'tGratingDirectionStepsPerOctave', ...
    'tLaserPeakStepsPerOctave', ...
    'tBaseGratingContrast', ...
    'tBaseGratingDirectionDeg', ...
    'tGratingContrast', ...
    'tGratingDirectionDeg', ...
    'tGratingElevationDeg', ...
    'tGratingAzimuthDeg', ...
    'tGratingHeightDeg', ...
    'tGratingWidthDeg', ...
    'tGratingSpatialFreqCPD', ...
    'tGratingSpeedDPS', ...
    'tGratingDurationMs', ...
    'tTrialLaserPowerMw', ...
    'tTrialLaserOnTimeMs', ...
    'tTrialLaserOffTimeMs', ...
    'tBlock2TrialNumber', ...
    'tNStimAccepted', ...
    'tLeverPressTimeMs', ...
    'tLeverReleaseTimeMs', ...
    'tCyclesOn', ...
    'nCyclesOn', ...
    'tTimeoutMs', ...
    'tStimOnMs', ...
    'targetStimOnMs', ...
    'tStimOffMs', ...
    'tTrialRemainingOffMs', ...
    'tTrialRemainingOnMs', ...
    'tTrialStartOn', ...
    'tTrialStartOff', ...
    'startReactState',...
    'stimTag',...
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

%calculate cycle duration
totalCycleTimeMs = input.stimOnTimeMs + input.stimOffTimeMs;
numberCyclesOn = input.nCyclesOn{trN};

%define leverpressoffset from wait
if input.tTrialStartOn{trN} == 1
    leverpressoffset = -(input.stimOnTimeMs-input.tTrialRemainingOnMs{trN}); 
elseif input.tTrialStartOff{trN} == 1
    leverpressoffset = input.tTrialRemainingOffMs{trN};
end
reqHoldTimeMs = double((totalCycleTimeMs*numberCyclesOn)+leverpressoffset);

input.leverpressoffset{trN} = leverpressoffset;

holdTimeMs = double((leverUpUs - leverDownUs)) / 1000;
reactTimeMs = (holdTimeMs - reqHoldTimeMs);

% add to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;
input.tTotalReqHoldTimeMs{trN} = reqHoldTimeMs;


% backward compat
input.gratingContrast = input.tGratingContrast;
input.laserPowerMw = input.tLaserPowerMw;
input.gratingDirectionDeg = input.tGratingDirectionDeg;

%Andrew's Post-Hoc Reaction Time Method
if input.targetStimOnMs{trN} <= 1;
    tCyclesShort = input.nCyclesOn{trN} - input.tCyclesOn{trN};
    input.targetStimOnMs{trN} = input.tStimOnMs{trN} + (totalCycleTimeMs*tCyclesShort);
end
input.postHocReactMs{trN} = input.tLeverReleaseTimeMs{trN} - input.targetStimOnMs{trN};


%% disp status
juiceAmtsMs = input.juiceTimesMsCell{trN};
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
stimStr = '';
if input.doVisualStim
    if input.doContrastDetect
        stimStr = [stimStr sprintf('ctst chg %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];
    elseif input.doOriDetect
        stimStr = [stimStr sprintf('direction %g%deg ', chop(input.gratingDirectionDeg{trN},2))];
    end
end
if input.doLaserStim
  stimStr = [stimStr sprintf('power %gmW ', ...
                    chop(input.laserPowerMw{trN}, 2))];
end
if ~input.doBlock2
  block2Str = '';
else
  block2Str = sprintf('b2tr %d ', input.tBlock2TrialNumber{trN});
end
itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
fprintf(1,'Hold %d, req %d, react %d, %s%s %s- %d rew %dms\n', ...
        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
        itiStr, ...
        stimStr, block2Str, ...
        nJ, round(juiceD));

%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotFlashingStim });



%% save variables for next trial
retval = input;

return


