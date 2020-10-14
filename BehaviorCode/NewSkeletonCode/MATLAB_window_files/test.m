function [retval] = MotorizedWheel(data_struct, input)
% Online Matlab processing for motorizedWheel.xml
% Code for Hull Lab - Duke Neurobiology
% Completed Oct 8th 2020 by Shuyang Jin
% Modified from MotorizedWheel
% call: exptSetupBridge, exptProcessBridgeInput, package data, then call subfunctions

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end
%Update for motorized wheel
ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/MatlabSharedCode');
oneValEachTrialNames = { ...
    'tTrialsDoneSinceStart', ...
    'tTrialStartMWTimestampMs', ...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'commandVoltage', ...
    'speed_mode', ...
    'voltage'};

exptSetupBridge;

 [input, eventsConsts, eventsTrial] ...
     = exptProcessBridgeInput(data_struct, input, oneValEachTrialNames);

trN = input.trialSinceReset;


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

