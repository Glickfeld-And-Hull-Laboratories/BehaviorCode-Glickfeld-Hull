function [retval] = MotorizedWheel(data_struct, input)
% Online Matlab processing for motorizedWheel.xml
% Code for Hull Lab - Duke Neurobiology
% Completed January 13th 2020 by Shuyang Jin
% Modified from CerebellarStim
% call: exptSetupBridge, exptProcessBridgeInput, package data, then call subfunctions

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end
%Update for motorized wheel
ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
oneValEachTrialNames = { ...
    'tCounter', ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'commandVoltage', ...
    'speed_mode', ...
    'voltage', ...
    'cTactileStimTurnedOn',...
    'cAuditoryStimTurnedOn',...
    'firstStimDelayNFrames',...
    'secondStimDelayNFrames',...
    'endDelayFrames'};

exptSetupBridge;

 [input, eventsConsts, eventsTrial] ...
     = exptProcessBridgeInput(data_struct, input, oneValEachTrialNames);

trN = input.trialSinceReset;

try
    input.counterTimesUs{trN} = mwGetEventTime(ds.events, ds.event_codec, 'counter', 'all', [], 1);
    input.counterValues{trN} = mwGetEventValue(ds.events, ds.event_codec, 'counter', 'all', 1) ;
catch
    input.counterTimesUs{trN} = NaN;
    input.counterValues{trN} = NaN;
end

try
    input.cAirpuffsOn{trN} = mwGetEventValue(ds.events, ds.event_codec, 'cTactileStimTurnedOn', 'all', 1) ;
catch
    input.cAirpuffsOn{trN} = NaN;
end

try
    input.cTonesOn{trN} = mwGetEventValue(ds.events, ds.event_codec, 'cAuditoryStimTurnedOn', 'all', 1) ;
catch
    input.cTonesOn{trN} = NaN;
end

%

%input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
%                              '~/Documents/MWorks/Data', ...
%                              mat2str(input.subjectNum), input.saveTime);
%save(input.savedDataName, 'input');
%  trStr = sprintf('Trial %1.0f Completed', trN);
% disp(trStr)

%% Run Plotting Function
% input.holdStartsMs = input.tThisTrialStartTimeMs;
input = exptRunSubfunctions(ds, input, { @plotCerebellarStim });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Returns calculated and passed variables for the next trial.
retval = input;
return

