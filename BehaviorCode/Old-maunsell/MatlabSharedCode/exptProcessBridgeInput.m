function [input eventsConsts eventsTrial] = exptProcessBridgeInput(data_struct, input, oneValEachTrialNames)
% This function processes event data and saves variables in input.
%
% Should be common to all experiments
% 
% MH 130105

ds = data_struct;
events = ds.events;
codecStruct = struct_vect2singleton(ds.event_codec);
ecs = exptConstants;


%% initial error checking
if ~exist('mwGetEventValue')
  error('Missing mwGetEventValue - is path correct?');
end


%% debug inputs
%dumpEvents(ds); 
%length(ds.events)
%save('/tmp/test.mat', 'data_struct', 'input');


%% setup vars on first call after pressing "play" in client
if nargin == 1 || ~isfield(input, 'trialSinceReset')
    % init input
    disp('First trial, initializing input');
    input.trialSinceReset = 1;
    input.startDateVec = datevec(now);
    input.savedEvents = {};
    input.eventCodecs = {};
    input.trialOutcomeCell = {};

    % init each vector with a blank cell
    nOne = length(oneValEachTrialNames);
    for iV = 1:nOne
        input.(oneValEachTrialNames{iV}) = {};
    end
else
  input.trialSinceReset = input.trialSinceReset+1;
end
thisTrialN = input.trialSinceReset;

%% misc check to be removed later
%assert(all(input.tooFastTimeMs >= 1));


%% save codec
if thisTrialN == 1
  input.eventCodecs{1} = ds.event_codec;
else
  lastCodecN = find(cellfun(@(x) ~isempty(x), input.eventCodecs), ...
                    1, 'last');
  if ~isequalwithequalnans(input.eventCodecs{lastCodecN}, ds.event_codec)
    error('bug: Codec seems to have changed - stopping!');
    %input.eventCodecs{thisTrialN} = ds.event_codec;
  end
end


%% process #announceStimulus events

asCodeN = codec_tag2code(input.eventCodecs{1}, '#announceStimulus');
codeV = cat(2, ds.events.event_code);
timeV = cat(2, ds.events.time_us);
asIx = codeV == asCodeN;
asTimeV = timeV(asIx);
%ds.events(asIx).data;

allData = {ds.events(asIx).data};
asActionV = cellfun(@(x) x.action, allData, 'UniformOutput', false);
asNameV = cellfun(@(x) x.name, allData, 'UniformOutput', false);
assert(all(strcmp(asActionV, 'draw')), ...
       'Found an action that was not ''draw'' - fix!');
unqNames = unique(asNameV);
for iN = 1:length(unqNames)
  tName = unqNames{iN};
  tNameIx = strcmp(asNameV, tName);
  tTimes = asTimeV(tNameIx);
  input.announceStimulusTimes(1).(tName){thisTrialN} = tTimes;
end
% remove  from savedEvents
events = events(~asIx);

% display some data about it
% $$$ at0 = input.announceStimulusTimes;
% $$$ wt0 = whos('at0');
% $$$ at1 = input.savedEvents;
% $$$ wt1 = whos('at1');
% $$$ fprintf(1, '--------------- savedEvents %.2f MB, asTimes %.2f MB\n', ...
% $$$         wt1.bytes/1e6, wt0.bytes/1e6);



%% save event stream
input.savedEvents{thisTrialN} = events;


%% are constants here?  if so, extract
nSync = length(mwGetEventValue(events, ds.event_codec, ...
                             'sync', 'all', 'ignoreMissing'));
eventsConsts = [];
allEvents = events;
if nSync >= 4 % is Justin Timberlake missing?
  % constants arrived.  Split out
  disp('Found dump of constants, splitting out from events');
  
  syncCode = codec_tag2code(ds.event_codec, 'sync');
  codeList = [events.event_code];
  syncNs = find(codeList == syncCode);

  % check a lot of assumptions.
  %   What should happen is sync=1, constants dumped included sync
  %   = 1, first trial starts, is interrupted by pressing stop,
  %   restart sends new codes with sync = 1, then sync = 0 stops
  %   trial.  Don't ask me why this is, pressing it once doesn't
  %   send the constants.
  errFlag = false;
  if ~ all(diff([events.time_us]) >= 0)
    disp('increasing in time?');     errFlag = true;
  elseif ~(syncNs(1) == 1)
    disp('sync structure'); errFlag = true;
  elseif ~(syncNs(end) == length(events))
    disp('last code'); errFlag = true;
  elseif ~(events(syncNs(1)).data == 1), 
    disp('err1'); errFlag = true; 
  elseif ~(events(syncNs(2)).data == 1), 
    disp('err2'); errFlag = true;
  elseif ~(events(syncNs(end-1)).data == 1), 
    disp('err3'); errFlag = true;
  elseif ~(events(syncNs(end)).data == 0), 
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

  % Take consts from first sync codes seen and use only last sync code run as start of trial
  % (Somewhat dangerous because we assume constants sent in that first trial)
  eventsConsts = events(1:syncNs(3));
  eventsTrial = events(syncNs(end-1):end);
  
  disp(sprintf('** Constants arriving in codec this trial! (tr %d)', ...
               input.trialSinceReset));
else
  eventsTrial = events;
end


%% make sure we're getting constants on first trial
if thisTrialN == 1

  % check that required variable names are present
  blank0 = codec_tag2code(ds.event_codec, 'subjectNum');  % will raise an error if missing
  blank0 = codec_tag2code(ds.event_codec, 'experimentXmlTrialId');  
    
  errFlag = false;
  if isempty(eventsConsts)
    errFlag = true;
  else
    tSubjectNumTemp = mwGetEventValue(eventsConsts, ds.event_codec, 'subjectNum', ...
                                      [], 'ignoremissing');
    if isempty(tSubjectNumTemp)
      errFlag = true;
    end
  end
  if errFlag
    error('** Error: must press start button twice on first trial!');
  end
end


%% process constants

trN = input.trialSinceReset;

% first find constList
persistentNames = codecStruct.tagname(codecStruct.persistant==1);
input.constList = persistentNames;
nConsts = length(input.constList);

%testing - verify against older code
%tcs = tempHADC8ConstantList;
%setdiff(persistentNames, tcs.constListDP)
%setdiff(tcs.constListDP, persistentNames)

% look for changes, all trials and load into constChangedonTrial
b={};
fieldList = fieldnames(input);
if input.trialSinceReset > 1
  for iC = 1:nConsts
    tCN = input.constList{iC};
    tV = mwGetEventValue(allEvents, ds.event_codec, tCN, 'last', 'ignoremissing');
    
    if ~isempty(tV)
      if ~any(ismember(fieldList, tCN)) || ~isequalwithequalnans(tV, input.(tCN))
        % different than last trial (or not a field - which only happens when adding vars to xml)
        % stuff into array saying what changed
        b(end+1,1:2) = {tCN, tV};
      end
    end
  end
end
input.constChangedOnTrial{input.trialSinceReset} = b;


% get constants from dump, first trial only  (eventsConsts)
for iC = 1:nConsts
  tCN = input.constList{iC};
  if input.trialSinceReset==1  % use dump on first trial
    tCV1 = mwGetEventValue(eventsConsts, ds.event_codec, tCN, [], ...
                          'ignoremissing');

    if isempty(tCV1)
      error(sprintf('*** Constant %s not found in dump - is it checked in Events box?!!!', ...
                   tCN));
    end
    input.(tCN) = tCV1;
  
  else  % if not trial 1, always copy const values in
    tCV = mwGetEventValue(events, ds.event_codec, tCN, 'last', ...
                          'ignoremissing');
    if ~isempty(tCV)
      input.(tCN) = tCV;
    end
  end
end

% make some into arrays
if any(ismember(input.constList, 'trPer80Level1'))
    input.trPer80V = [input.trPer80Level1 input.trPer80Level2 input.trPer80Level3 ...
                      input.trPer80Level4 input.trPer80Level5 input.trPer80Level6 ...
                      input.trPer80Level7 input.trPer80Level8];
    input.trPer80History{trN} = input.trPer80V;
end
if any(ismember(input.constList, 'block2TrPer80Level1'))
    input.block2TrPer80V = [input.block2TrPer80Level1 input.block2TrPer80Level2 input.block2TrPer80Level3 ...
                        input.block2TrPer80Level4 input.block2TrPer80Level5 input.block2TrPer80Level6 ...
                        input.block2TrPer80Level7 input.block2TrPer80Level8];              
    input.block2TrPer80History{trN} = input.block2TrPer80V;
end

%% process trial outcome
thisExptOutcomes = intersect(ecs.trialOutcomes, codecStruct.tagname);
nOutcomes = length(thisExptOutcomes);
nFoundThisTrial = repmat(NaN, [1 nOutcomes]);
for iO = 1:nOutcomes
    nFoundThisTrial(iO) ...
        = length(mwGetEventValue(eventsTrial, ds.event_codec, ...
                            thisExptOutcomes{iO}, 'all', 'ignoreMissing'));
end

if sum(nFoundThisTrial) ~= 1 
    %disp(eventsTrial); disp(ds.event_codec); disp(nSync);
    nFoundThisTrial
    error('error - in outcome calculation, check'); 
end
outcomeStr = thisExptOutcomes{nFoundThisTrial == 1};

input.trialOutcomeCell{trN} = outcomeStr;



%% save first trial status
if input.trialSinceReset == 1
  input.firstTrConsts = input;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process vectors, one value each trial

for iV = 1:length(oneValEachTrialNames)
    tName = oneValEachTrialNames{iV};
    input.(tName){trN} = mwGetEventValue(eventsTrial, ds.event_codec, tName, 'last', 1);
end

% total reward times
juiceAmtsUs = mwGetEventValue(eventsTrial, ds.event_codec, 'juice', 'all');
juiceAmtsMs = double(juiceAmtsUs(juiceAmtsUs~=0))/1000;
input.juiceTimesMsCell{thisTrialN} = juiceAmtsMs;




%% return outputs
retval = input;
return

