function [retval] = CerebellarStim(data_struct, input)
% Online Matlab processing for cerebellarStim.xml
% Code for Hull Lab - Duke Neurobiology
% Completed August 29th, 2014 by Andrew McKinney
% call: exptSetupBridge, exptProcessBridgeInput, package data, then call subfunctions

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
oneValEachTrialNames = { ...
    'tCounter', ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'cTactileStimTurnedOn', ...
    'tResistanceStartMS',...
    'tReverseVStimTimeMs',...
    'cReverse',...
    'tReverseDelayMS',...
    'tDotDirectionDeg',...
    };

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

