
    'laserTrigger', ...    'laserTrigger', ...function [retval] = HoldAndDetectConstant(data_struct, input)
% Main matlab online function for HADC8
%
%  MH 100115: created
%  MH 130107: refactored into new expt* functions.
% call: exptSetupBridge, exptProcessBridgeInput, do your own thing, then exptRunSubfunctions

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/Old-maunsell/MatlabSharedCode');

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
    'tGratingMaxSpeedStepDPS', ...
    'tGratingSpeedStepsPerOctave', ...
    'tLaserPeakStepsPerOctave', ...
    'tBaseGratingContrast', ...
    'tBaseGratingDirectionDeg', ...
    'tBaseGratingSpeedDPS', ...
    'tGratingContrast', ...
    'tGratingDirectionDeg', ...
    'tGratingElevationDeg', ...
    'tGratingAzimuthDeg', ...
    'tGratingHeightDeg', ...
    'tGratingWidthDeg', ...
    'tGratingSpatialFreqCPD', ...
    'tGratingSpeedDPS', ...
    'tGratingStartingPhaseDeg', ...
    'tGratingDurationMs', ...
    'tTrialLaserPowerMw', ...
    'tTrialLaserOnTimeMs', ...
    'tTrialLaserOffTimeMs', ...
    'tDoAllCorrects', ...
    'tBlock2TrialNumber', ...
    'tDoNoStimulusChange', ...
    'tNStimAccepted', ...
    'laserTrigger', ...
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


% backward compat
input.gratingContrast = input.tGratingContrast;
input.laserPowerMw = input.tLaserPowerMw;
input.gratingDirectionDeg = input.tGratingDirectionDeg;


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
    elseif input.doSpeedDetect
        stimStr = [stimStr sprintf('speed %g%dps ', chop(abs((double(input.tGratingSpeedDPS{trN}) - double(input.tBaseGratingSpeedDPS{trN}))),2))];
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
input = exptRunSubfunctions(ds, input, { @plotOnlineHist });



%% save variables for next trial
retval = input;

return


