function [retval] = HoldAndDetectBeta(data_struct, input)

% This function processes event data and saves variables in input.
%
% Then it calls plot functions in a
% try-catch block so that errors in plot functions don't affect
% variable saving.
%  MH 100115
% $Id$

beep off;  % otherwise played through speakers to animal!
format compact; 

ds = data_struct;

% debug
dumpEvents(ds); 
%length(ds.events)
%save('/tmp/test.mat', 'data_struct', 'input');

% First call after pressing "play" in client
if nargin == 1 || ~isfield(input, 'trialSinceReset')
    % init input
    disp('First trial, initializing input');
    input.trialSinceReset = 1;
    input.reactTimesMs = {};
    input.holdTimesMs = {};
    input.reqHoldTimeMs = {};
    input.tooFastTimeMs = [];
    input.trialOutcomeCell = {};
    input.holdStartsMs = {};
    input.juiceTimesMsCell = {};
    input.constChangedOnTrial = {};
else
  input.trialSinceReset = input.trialSinceReset+1;
  assert(all(input.tooFastTimeMs >= 1));
end
thisTrialN = input.trialSinceReset;

%% quick path check 
if ~exist('mwGetEventValue')
  error('Missing mwGetEventValue - is path correct?');
end


%% are constants here?  if so, extract
nSync = length(mwGetEventValue(ds.events, ds.event_codec, ...
                             'sync', 'all', 'ignoreMissing'));
constsStruct = [];
if nSync >= 4 % is Justin Timberlake missing?
  % constants arrived.  Split out
  disp('Found dump of constants, splitting out from events');
  
  syncCode = codec_tag2code(ds.event_codec, 'sync');
  codeList = [ds.events.event_code];
  syncNs = find(codeList == syncCode);

  % check a lot of assumptions.
  %   What should happen is sync=1, constants dumped included sync
  %   = 1, first trial starts, is interrupted by pressing stop,
  %   restart sends new codes with sync = 1, then sync = 0 stops
  %   trial.  Don't ask me why this is, pressing it once doesn't
  %   send the constants.
  errFlag = false;
  if ~ all(diff([ds.events.time_us]) >= 0)
    disp('increasing in time?');     errFlag = true;
  elseif ~(syncNs(1) == 1)
    disp('sync structure'); errFlag = true;
  elseif ~(syncNs(end) == length(ds.events))
    disp('last code'); errFlag = true;
  elseif ~(ds.events(syncNs(1)).data == 1), 
    disp('err1'); errFlag = true; 
  elseif ~(ds.events(syncNs(2)).data == 1), 
    disp('err2'); errFlag = true;
  elseif ~(ds.events(syncNs(end-1)).data == 1), 
    disp('err3'); errFlag = true;
  elseif ~(ds.events(syncNs(end)).data == 0), 
    disp('err4'); errFlag = true;
  end

  if errFlag == true
    disp('** error in sync code pattern: check and fix code');
    disp(syncNs);
    disp(nSync);
    retval = []
    return
  end
  
  if nSync > 4
    disp('Too many sync codes:');
    disp(nSync);
    % just display a message
    %   errFlag = true;
  end

  % if 6 codes, pull out consts from first block and leave last
  % (Somewhat dangerous)
  constsStruct = ds.events(1:syncNs(3));
  trialStruct = ds.events(syncNs(end-1):end);
  ds.events = trialStruct;
end


%% process constants
input.constList = { 'tooFastTimeMs', 'subjectNum', 'doWaitForUp', ...
                    'doLeverSolenoidAllTrials', ...
                    'doHoldTone', 'reactTimeMs', 'randReqHoldMaxMs', ...
                    'doExtendItiOnShortPrevTrial', ...
                    'fixedReqHoldTimeMs', 'earlyTimeoutMs', ...
                    'missedTimeoutMs', ...
                    'itiTimeMs', 'postRewardMs', ...
                    'minRewardUs', 'maxRewardUs', ...
                    'maxConsecCorrects', 'jackpotProb', ...
                    'jackpotRewardSizeUs', 'doLever' };
nConsts = length(input.constList);
% get constants from dump, first trial only  (constsStruct)

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



% look for changes, all trials  (ds.events)
b={};
for iC = 1:nConsts
  tCN = input.constList{iC};
  tV = mwGetEventValue(ds.events, ds.event_codec, tCN, [], 'ignoremissing');

  if ~isempty(tV), 
    % stuff into array saying what changed
    b(end+1,1:2) = {tCN, tV};
  end
end
input.constChangedOnTrial{input.trialSinceReset} = b;


%% make sure we're getting constants on first trial
if ~isfield(input, 'subjectNum') || isempty(input.subjectNum)

  disp('** Error: must press start button twice on first trial!');
  retval = []; 
  return
end
if ~isempty(mwGetEventValue(ds.events, ds.event_codec, 'subjectNum', ...
                            [], 'ignoremissing'))
  disp(sprintf('** Constants arriving in codec this trial! (tr %d)', ...
               input.trialSinceReset));
end

%% process corrects
nS = length(mwGetEventValue(ds.events, ds.event_codec, ...
                             'success', 'all', 'ignoreMissing'));
nF = length(mwGetEventValue(ds.events, ds.event_codec, ...
                             'failure', 'all', 'ignoreMissing'));
nI = length(mwGetEventValue(ds.events, ds.event_codec, ...
                             'ignore', 'all', 'ignoreMissing'));
if sum(nS+nF+nI) > 1, 
  disp(ds.events); disp(ds.event_codec); disp(nSync);
  disp(nS); disp(nF); disp(nI);
  error('error - 1, check'); 
end
if nS == 1, tStr = 'success'; 
elseif nF == 1, tStr = 'failure';
elseif nI == 1, tStr = 'ignore';
else
  disp('Error!  Missing trial outcome variable this trial');
  tStr = 'error-missing';
end
% save
input.trialOutcomeCell{input.trialSinceReset} = tStr;

%% process reaction times for this trial
codes = [ds.events.event_code];

stimOnUs = mwGetEventTime(ds.events, ds.event_codec, 'stimulusOn', 1);
totalHoldTimeMs = mwGetEventValue(ds.events, ds.event_codec, 'tTotalReqHoldTimeMs');
%trialStartUs = mwGetEventTime(ds.events, ds.event_codec, 'trialStart', 1);
leverDownUs = mwGetEventTime(ds.events, ds.event_codec, 'leverResult', 1, 1);
leverUpUs = mwGetEventTime(ds.events, ds.event_codec, 'leverResult', 2, 0);

%reactTimeMs = (leverUpUs - stimOnUs) / 1000;
holdTimeMs = (leverUpUs - leverDownUs) / 1000;
reactTimeMs = (holdTimeMs - totalHoldTimeMs);

% add to array
input.holdStartsMs{input.trialSinceReset} = leverDownUs/1000;
input.holdTimesMs{input.trialSinceReset} = holdTimeMs;
input.reactTimesMs{input.trialSinceReset} = reactTimeMs;
input.reqHoldTimeMs{input.trialSinceReset} = totalHoldTimeMs;

% total reward times
juiceAmtsUs = mwGetEventValue(ds.events, ds.event_codec, 'juice', 'all');
juiceAmtsMs = juiceAmtsUs(juiceAmtsUs~=0)/1000;
input.juiceTimesMsCell{thisTrialN} = juiceAmtsMs;

% disp status
juiceD = min(juiceAmtsMs);  % if isempty, this will be empty
nJ = length(juiceAmtsMs);
fprintf(1, 'Hold %4d, req %4d, react %4d, reward (%1d) %3d\n ', ...
        round(holdTimeMs), round(totalHoldTimeMs), round(reactTimeMs), nJ, round(juiceD));
        
%% save variables for next trial
retval = input;

%% run subfunctions

addpath('/Library/Application Support/MWorks/Scripting/Matlab/tools-mh');

tic
try
  input = saveMatlabState(data_struct, input);

  % do the plot
  if ~isfield(input, 'lastPlotElapsedS')
    input.lastPlotElapsedS = 0;
  end
  if input.lastPlotElapsedS > 2
    % skip every other plot to make sure we catch up
    disp(sprintf('Skipping plot, last elapsed was %5.1fs', ...
                 input.lastPlotElapsedS));
        input.lastPlotElapsedS = 0;
  else
    % last one was fast, plot this time
    tic
    input = plotOnlineHist(ds, input);  % this gets the CHANGED event stream w/ first trial consts extracted
    input.lastPlotElapsedS = toc;
    %disp(input.lastPlotElapsedS);
  end
    

  %input = testUpload(data_struct, input);

catch ex
  disp('??? Error in subfunction; still saving variables for next trial')
  printErrorStack(ex);
end
toc
%% save variables for next trial
retval = input;

return

