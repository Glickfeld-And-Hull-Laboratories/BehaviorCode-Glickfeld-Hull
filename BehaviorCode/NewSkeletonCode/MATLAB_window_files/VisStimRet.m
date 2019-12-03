function [retval] = VisStimRet(data_struct, input)
%Main MATLAB interface function for VisStimRet.xml.
%   
if nargin < 2,
    input = [];
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');

% Vars set on each trial
oneValEachTrialNames = { ...
'tTempStim', ...
'counter', ...
'tTrialsDoneSinceStart' , ...
'tTrialStartMWTimestampMs', ...
'tThisTrialStartTimeMs' , ...
'tLastTrialStartTimeMs'   , ...  
'rrStimulusNumber' , ...
'tStimulusNumber', ...
'tBaseGratingContrast',...
'tBaseGratingDirectionDeg',...
'tGratingAzimuthDeg', ...
'tGratingElevationDeg', ...
'tGratingDirectionDeg', ...
'tGratingContrast' , ...
'tGratingDiameterDeg', ...
'tGratingSpatialFreqCPD', ...
'tGratingTemporalFreqCPS', ...
'tGratingStartingPhaseDeg', ...
'tGratingSpeedDPS', ...
'tGratingLoomSpeedDPS', ...
'tAnnulusGratingDiameterDeg',...
'tGreymaskDiameterDeg',...
'tDotSpeedDPS', ...
'tDotSizeDeg', ...
'tDotDirectionDeg',...
'tDotContrast',...
'tDotCoherence',...
'tDotDensity',...
'tDotElevationDeg',...
'tDotAzimuthDeg',...
'tDotFieldSizeDeg',...
'tNStimAccepted', ...
'frameCountN', ...
'tGratingTemporalFreqCPS', ...
'tTrialLaserPowerMw', ...
'tTrialLaserPowerMw_trigger',...
'setDotSpeedDPS'};


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
    'experimentXmlTrialId', ...
    'nScansOn', ...
    'nScansOff', ...
    'stopAfterNTrials', ...
    'doRand',...
    'doMatrix',...
    'doRetStim',...
    'doConStim',...
    'doSizeStim',...
    'doAnnulusStim',...
    'doAnnulusSizeStim',...
    'doSFStim',...
    'doTFStim',...
    'doDirStim',...
    'doLEDStim',...
    'doPhaseStim',...
    'doLoomStim',...
    'doGaussianMask',...
    'doEllipseMask',...
    'doRunFeedback',...
    'doMovingDots',...
    'doDotsSpeedStim',...
    'doDotsSizeStim',...
    'doDotsCoherenceStim',...
    'doDotsDirectionStim',...
    'doDotsContrastStim',...
    'doDotsDensityStim',...
    'runFeedbackGain',...
    'dotAzimuthDeg',...
    'dotElevationDeg',...
    'nTotalImagingFrames',...
    'frameImagingRateMs',...
    'frameImagingExposureMs',...
    'frameImagingBinning',...
    'baseGratingContrast',...
    'baseGratingDirectionDeg',...
    'gratingElevationDeg',...
    'gratingElevationStepDeg',...
    'gratingElevationStepN',...
    'gratingAzimuthDeg',...
    'gratingAzimuthStepDeg',...
    'gratingAzimuthStepN',...
    'gratingContrast',...
    'gratingContrastStepLog',...
    'gratingContrastStepDir',...
    'gratingContrastStepN',...
    'gratingDirectionDeg',...
    'gratingDirectionStepDeg',...
    'gratingDirectionStepN',...
    'gratingDiameterDeg',...
    'gratingDiameterStepDeg',...
    'gratingDiameterStepN',...
    'gratingSpatialFreqCPD',...
    'gratingSpatialFreqStepLog',...
    'gratingSpatialFreqStepDir',...
    'gratingSpatialFreqStepN',...
    'gratingTemporalFreqCPS',...
    'gratingTemporalFreqStepLog',...
    'gratingTemporalFreqStepDir',...
    'gratingTemporalFreqStepN',...
    'gratingStartingPhaseDeg',...
    'gratingStartingPhaseStepDeg',...
    'gratingStartingPhaseStepN',...
    'greymaskDiameterDeg',...
    'AnnulusgratingDiameterStepN',...
    'AnnulusgratingDiameterMinN',...
    'AnnulusgratingDiameterMaxN',...
    'loomSpeedDPS',...
    'maxLoomDiameterDeg',...
    'loomSpeedStepLog',...
    'loomSpeedStepDir',...
    'loomSpeedStepN',...
    'dotElevationDeg',...
    'dotAzimuthDeg',...
    'dotFieldSizeDeg',...
    'dotContrast',...
    'dotContrastStepLog',...
    'dotContrastStepDir',...
    'dotContrastStepN',...
    'dotCoherence',...
    'dotCoherenceStepLog',...
    'dotCoherenceStepDir',...
    'dotCoherenceStepN',...
    'dotDirectionDeg',...
    'dotDirectionStepDeg',...
    'dotDirectionStepN',...
    'dotSizeDeg',...
    'dotSizeStepDeg',...
    'dotSizeStepN',...
    'dotSpeedDPS',...
    'dotSpeedStepLog',...
    'dotSpeedStepDir',...
    'dotSpeedStepN',...
    'dotDensity',...
    'dotDensityStep',...
    'dotDensityStepN',...
    'doRand',...
    'doWheelSpeed'};


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

%% Counter/Frames Synchronization Section
try
    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

%% Running events Section
try
    input.wheelSpeedTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'wheelSpeed', 'all', [], 1);
    input.wheelSpeedValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'wheelSpeed', 'all', 1) ;
catch
    input.wheelSpeedTimesUs{trN} = NaN;
    input.wheelSpeedValues{trN} = NaN;
end

%% Speed changes if feedback
if input.doRunFeedback
  try
      input.dotSpeedTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'setDotSpeedDPS', 'all', [], 1);
      input.dotSpeedValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'setDotSpeedDPS', 'all', 1) ;
  catch
      input.dotSpeedTimesUs{trN} = NaN;
      input.dotSpeedValues{trN} = NaN;
  end
end


input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');

if ~input.doMovingDots
disp(sprintf('Trial %d: Az %d, El %d, Con %1.2f, Dir %d deg, SF %1.2f CPD, TF %d Hz, Sz %d deg', ...
    trN, input.tGratingAzimuthDeg{trN}, input.tGratingElevationDeg{trN}, input.tGratingContrast{trN}, input.tGratingDirectionDeg{trN}, input.tGratingSpatialFreqCPD{trN}, input.tGratingTemporalFreqCPS{trN}, chop(input.tGratingDiameterDeg{trN},2)));
elseif input.doMovingDots & ~input.doRunFeedback
  disp(sprintf('Trial %d: Az %d, El %d, Con %1.2f, Direction %d deg, Speed %d dps, Coh %1.2f', ...
    trN, input.dotAzimuthDeg(1), input.dotElevationDeg(1), input.tDotContrast{trN}, input.tDotDirectionDeg{trN}, input.tDotSpeedDPS{trN}, input.tDotCoherence{trN}));
elseif input.doMovingDots & input.doRunFeedback
  disp(sprintf('Trial %d: Az %d, El %d, Con %1.2f, Direction %d deg, Speed: feedback, Coh %1.2f', ...
    trN, input.tDotAzimuthDeg(1), input.tDotElevationDeg(1), input.tDotContrast{trN}, input.tDotDirectionDeg{trN}, input.tDotCoherence{trN}));
end

%% Run Plotting Function
input.holdStartsMs = input.tThisTrialStartTimeMs;
input = exptRunSubfunctions(ds, input, { @plotVisStim });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%trN = input.trialSinceReset;

%stimStr = [stimStr sprintf('Contrast %g%% ', chop(abs(double(input.gratingContrast{trN})-double(input.tBaseGratingContrast{trN})),2)*100)];

%% Returns calculated and passed variables for the next trial.
retval = input;
return

