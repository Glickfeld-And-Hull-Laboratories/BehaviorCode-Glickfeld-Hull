function [retval] = HoldAndDetect_2P_Frames(data_struct, input)
% Main matlab online function for HoldAndDetect_Frames.xml
%
% Adapted by Andrew McKinney, finalized January 20th, 2015

if nargin < 2,
    input = [];
    input.saveTime = datestr(now, 'yymmdd-HHMM'); 
end

ds = data_struct;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/Old-maunsell/MatlabSharedCode');

varsOneValEachTrial = { ...
    'tItiTimeMs',...
    'tThisTrialStartTimeMs', ...
    'tLastTrialStartTimeMs', ...
    'tTrialStartMWTimestampMs', ...
    'tRewardTrial', ...
    'tGratingDirectionDeg', ...
    'dGratingDirectionDeg', ...
};

exptSetupBridge;

[input eventsConsts eventsTrial] ...
    = exptProcessBridgeInput(data_struct, input, varsOneValEachTrial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trN = input.trialSinceReset;

%% disp status
juiceAmtsMs = input.juiceTimesMsCell{trN};
juiceD = sum(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
stimStr = '';
stimStr = [stimStr sprintf('direction %g%deg ', chop(input.gratingDirectionDeg{trN},2))];

fprintf(1,%s- rew %dms\n', ...
        stimStr, round(juiceD));


%% Counter/Frames Synchronization Section
%try
%    input.counterTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'counter', 'all', [], 1);
%    input.counterValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'counter', 'all', 1) ;
%catch
%    input.counterTimesUs{trN} = NaN;
%    input.counterValues{trN} = NaN;
%end
%
%try
%    input.quadratureTimesUs{trN} = mwGetEventTime(eventsTrial, ds.event_codec, 'quadrature', 'all', [], 1);
%    input.quadratureValues{trN} = mwGetEventValue(eventsTrial, ds.event_codec, 'quadrature', 'all', 1) ;
%catch
%    input.quadratureTimesUs{trN} = NaN;
%    input.quadratureValues{trN} = NaN;
%end

%% Mati's Lever Tracking Addition -- Polls for all lever changes (high or low)

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


%% run subfunctions
input = exptRunSubfunctions(ds, input, { @plotOnlineHist });



%% save variables for next trial
retval = input;

return


