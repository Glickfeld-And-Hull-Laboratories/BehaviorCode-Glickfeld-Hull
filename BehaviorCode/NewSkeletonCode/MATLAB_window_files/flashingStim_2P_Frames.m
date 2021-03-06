function [retval] = HoldAndDetectConstant(data_struct, input)
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
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');

varsOneValEachTrial = {...
    'tTotalRewardTimeUs', ...
    'tFakeMouseReactFrames', ...
    'tStartTrialWaitForPressTimeMs', ...
    'tConsecErrors', ...
    'tConsecTimeoutStartTime', ...
    'tNRewards', ...
    'tInterRewardIntervalMs', ...
    'tRewardAddPerMsHoldUs', ...
    'tStimTurnedOn', ...
    'tItiWaitTimeMs', ...
    'tItiWaitFrames', ...
    'tTotalReqHoldFrames', ...
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
    'tTrialLaserPowerMw', ...
    'tTrialLaserOnTimeMs', ...
    'tTrialLaserOffTimeMs', ...
    'tBlock2TrialNumber', ...
    'tNStimAccepted', ...
    'tLeverPressTimeMs', ...
    'tLeverReleaseTimeMs', ...
    'tDoAuditoryDetect', ...
    'tCyclesOn', ...
    'nCyclesOn', ...
    'tStimOnMs', ...
    'tStimTurnedOn',...
    'targetStimOnMs', ...
    'tStimOffMs', ...
    'startReactState',...
    'stimTag',...
    'isTooFast',...
    'isCatchTrial',...
    'cItiStart',...
    'cLeverDown',...
    'cFirstStim',...
    'cStimOn',...
    'cStimOff',...
    'cTargetOn',...
    'cEndOn',...
    'cEndOff',...
    'cLeverUp',...
    'cAuditoryStimOn',...
    'catchCyclesOn',...
    'tShortCatchTrial',...
    'tFalseAlarm',...
    'cCatchOn',...
    'tCatchGratingDirectionDeg',...
    'tCatchGratingContrast',...
    'tSoundTargetAmplitude',...
    'tSoundCatchAmplitude',...
    'tStimOnTimeMs',...
    'tStimOffTimeMs',...
    'nFramesOff',...
    'nFramesOn',...
    'tTotalStimFrames',...
    'multVal',...
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
holdTimeMs = double((leverUpUs - leverDownUs)) / 1000;

if isfield(input, 'doRandStimOffTime') 
    if input.doRandStimOffTime == 1
        tempTimes = mwGetEventValue(eventsTrial, ds.event_codec, 'nFramesOff', 'all', 1);
        input.tFramesOff{trN} = tempTimes(2:end);
        totalStimTimeMs= (input.tTotalStimFrames{trN}./double(input.frameRateHz)).*1000;
        numberCyclesOn = input.nCyclesOn{trN};
        nCyclesRemaining = numberCyclesOn - input.tCyclesOn{trN};
        avgCycleTime = ((input.stimOffTimeMs+input.minStimOffTimeMs)./2) + input.stimOnTimeMs;
        reqHoldTimeMs = totalStimTimeMs + (nCyclesRemaining*avgCycleTime);
        reactTimeMs = double(holdTimeMs - reqHoldTimeMs);
    else
        totalCycleTimeMs = ((input.nFramesOn{1} + input.nFramesOff{1})./double(input.frameRateHz)).*1000;
        numberCyclesOn = input.nCyclesOn{trN};
        reqHoldTimeMs = double(totalCycleTimeMs*numberCyclesOn);
        holdTimeMs = double((leverUpUs - leverDownUs)) / 1000;
        reactTimeMs = (holdTimeMs - reqHoldTimeMs);
    end
else
    totalCycleTimeMs = ((input.nFramesOn + input.nFramesOff)./double(input.frameRateHz)).*1000;
    numberCyclesOn = input.nCyclesOn{trN};
    reqHoldTimeMs = double(totalCycleTimeMs*numberCyclesOn);
    reactTimeMs = (holdTimeMs - reqHoldTimeMs);
end

%find StimOn frames
stimOnFrames = zeros(1, input.maxCyclesOn);
 for istim = 1:numberCyclesOn
    try
        stimOnFrames(1,istim) = mwGetEventValue(eventsTrial, ds.event_codec, 'cStimOn', istim);
    catch err
        stimOnFrames(1,istim) - NaN;
    end
 end

%track contrast by presentation
if input.doRandCon
    multVals = mwGetEventValue(eventsTrial, ds.event_codec, 'multVal', 'all', 1);
    input.tBaseGratingContrast{trN} = double(input.baseGratingContrast) + (double(multVals) .* double(input.conDiff));
end

auditoryStimOnFrames = zeros(1, input.maxCyclesOn);
for istim = 1:numberCyclesOn
    try
        auditoryStimOnFrames(1,istim) = mwGetEventValue(eventsTrial, ds.event_codec, 'cAuditoryStimOn', istim);auditoryStimOnFrames(1,istim) = mwGetEventValue(eventsTrial, ds.event_codec, 'cAuditoryStimOn', istim);
    catch err
        auditoryStimOnFrames(1,istim) = NaN;
    end
end

% add to array
input.holdStartsMs{trN} = leverDownUs/1000;
input.holdTimesMs{trN} = holdTimeMs;
input.reactTimesMs{trN} = reactTimeMs;
input.tTotalReqHoldTimeMs{trN} = reqHoldTimeMs;
input.cStimOnFrames{trN} = stimOnFrames;
input.cAuditoryStimOnFrames{trN} = auditoryStimOnFrames;


% backward compat
input.gratingContrast = input.tGratingContrast;
input.laserPowerMw = input.tLaserPowerMw;
input.gratingDirectionDeg = input.tGratingDirectionDeg;

% short catch trial calculations

  if input.tShortCatchTrial{trN}
    if input.tFalseAlarm{trN}
        input.catchTrialOutcomeCell{trN} = 'FA';
    end
    if isempty(input.cCatchOn{trN})
        input.cCatchOn{trN} = NaN;
        input.catchTrialOutcomeCell{trN} = 'failure';
    end
    if (input.cLeverUp{trN}-input.cCatchOn{trN})>input.nFramesReact
        input.catchTrialOutcomeCell{trN} = 'CR';
    end
    if (input.cLeverUp{trN}-input.cCatchOn{trN})<input.nFramesTooFast
        input.catchTrialOutcomeCell{trN} = 'failure';
    else
        input.catchTrialOutcomeCell{trN} = 'NaN';
    end
  end



%% Counter/Frames Synchronization Section
try
    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

try
    input.quadratureTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'quadrature', 'all', [], 1);
    input.quadratureValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'quadrature', 'all', 1) ;
catch
    input.quadratureTimesUs{trN} = NaN;
    input.quadratureValues{trN} = NaN;
end

try
    input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
    input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
catch
    input.wheelSpeedTimesUs{trN} = NaN;
    input.wheelSpeedValues{trN} = NaN;
end

%% disp status
juiceAmtsMs = input.juiceTimesMsCell{trN};
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
stimStr = '';
if input.tBlock2TrialNumber{trN} == 0
    if input.doVisualStim
        if input.doContrastDetect
            stimStr = [stimStr sprintf('ctst chg %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];
        elseif input.doOriDetect
            stimStr = [stimStr sprintf('direction %g%deg ', chop(input.gratingDirectionDeg{trN},2))];
        end
    end
    if input.doAuditoryStim
        if input.doAuditoryDetect
            % stimStr = [stimStr sprintf('pitch chg %g%MHz ', chop(double(input.tTonePitchMHz{trN})- double(input.baseTonePitchMHz{trN}),2))];
        end
    end
    if input.doLaserStim
      stimStr = [stimStr sprintf('power %gmW ', ...
                        chop(input.laserPowerMw{trN}, 2))];
    end
else
    if input.block2DoVisualStim
        if input.block2DoContrastDetect
            stimStr = [stimStr sprintf('ctst chg %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];
        elseif input.block2DoOriDetect
            stimStr = [stimStr sprintf('direction %g%deg ', chop(input.gratingDirectionDeg{trN},2))];
        elseif input.doContrastDetect
            stimStr = [stimStr sprintf('ctst chg %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];
        elseif input.doOriDetect
            stimStr = [stimStr sprintf('direction %g%deg ', chop(input.gratingDirectionDeg{trN},2))];
        end
    end
    if input.block2DoAuditoryStim
        if input.block2DoAuditoryDetect
            stimStr = [stimStr sprintf(' tone chg')];
            % future compat for tone steps
            % stimStr = [stimStr sprintf('pitch chg %g%MHz ', chop(double(input.tTonePitchMHz{trN})- double(input.tBaseTonePitchMHz{trN}),2))];
        end
    end
end

if ~input.doBlock2
  block2Str = '';
else
  block2Str = sprintf(' b2tr %d ', input.tBlock2TrialNumber{trN});
end
itiStr = sprintf('iti %d, ', round(input.tItiWaitTimeMs{trN}));
fprintf(1,'Hold %d, req %d, react %d, %s %s- %d rew %dms\n', ...
        round(holdTimeMs), round(input.tTotalReqHoldTimeMs{trN}), round(reactTimeMs), ...
        stimStr, block2Str, ...
        nJ, round(juiceD));

%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotFlashingStim });



%% save variables for next trial
retval = input;

return


