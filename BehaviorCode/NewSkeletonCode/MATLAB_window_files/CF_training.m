function [retval] = CF_training(data_struct, input)
% Main matlab online function for CF_training
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
    'tThisTrialStartTimeMs', ...
};

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

input.savedDataName = sprintf('%s/data-i%03s-%s.mat', ...
                              '~/Documents/MWorks/Data', ...
                              mat2str(input.subjectNum), input.saveTime);
save(input.savedDataName, 'input');

%% save variables for next trial
retval = input;

return


