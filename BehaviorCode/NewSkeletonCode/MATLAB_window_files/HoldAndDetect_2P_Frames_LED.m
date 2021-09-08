function [retval] = HoldAndDetect_2P_Frames(data_struct, input)
% Main matlab online function for HoldAndDetect_Frames.xml
%
% Adapted by Andrew McKinney, finalized January 20th, 2015

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');

varsOneValEachTrial = { ...
    ...%'holdStartsMs', ...
    'tTotalReqHoldTimeMs', ...
    'tTotalRewardTimeUs', ...
    'tRandReqHoldTimeMs', ...
    'tTotalReqHoldFrames', ...
    'tItiWaitFrames',...
    'tFakeMouseReactMs', ...
    'tFakeMousePressMs', ...
    'tFakeMouseReactFrames', ...
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
    'tDoAllCorrects', ...
    'tBlock2TrialNumber', ...
    'tDoNoStimulusChange', ...
    'tNStimAccepted', ...
    'tRewardOmissionTrial',...
    'tRewardDelayTrial',...
    'tItiUnexpectedRewardTrial',...
    'tRewardDelayDurationMs',...
    'tDoNoStimulusChange',...
    'cTrialStart', ...
    'cTrialEnd', ...
    'cLeverDown', ...
    'cTargetOn', ...
    'cLeverUp', ...
    'soundTargetAmplitude',...
};

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% process reaction times for this trial

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


%% Counter/Frames Synchronization Section
try
    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

try
    input.quadratureTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'quadrature', 'all', [], 1);
    input.quadratureValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'quadrature', 'all', 1) ;
catch
    input.quadratureTimesUs{trN} = NaN;
    input.quadratureValues{trN} = NaN;
end

%% Mati's Lever Tracking Addition -- Polls for all lever changes (high or low)
try
    input.leverTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'FIO1', 'all', [], 1);
    input.leverValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'FIO1', 'all', 1) ;
catch
    input.leverTimesUs{trN} = NaN;
    input.leverValues{trN} = NaN;
end

%% LG adds to keep track of stimulus change time
try
    input.stimOnUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'stimulusOn', 1);
catch
    input.stimOnUs{trN} = NaN;
end

%% Adds lickometer tracking for counter2
try
    input.lickometerTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter2', 'all', [], 1);
    input.lickometerValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter2', 'all', 1) ;
catch
    input.lickometerTimesUs{trN} = NaN;
    input.lickometerValues{trN} = NaN;
end

% if input.doLED==1
%    input.LEDTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'laserTriggerFIO', 'all', [], 1);
% else
% input.LEDTimesUs{trN} = NaN;
% end



%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotOnlineHist });



%% save variables for next trial
retval = input;

return


