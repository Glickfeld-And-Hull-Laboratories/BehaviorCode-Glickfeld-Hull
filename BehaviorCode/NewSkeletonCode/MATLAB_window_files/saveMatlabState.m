function input = saveMatlabState(data_struct, input)

% essential consts
cs = holdanddetect_constants;

% compute whether to do export on this trial
nTrialsToSkip = 6;
maxSToSkip = 100;

%% initialize variables
nTrial = input.trialSinceReset;
thisTimeMs = input.holdStartsMs{nTrial};
if nTrial == 1 || ~isfield(input, 'lastTimeSavedMs')
  input.lastTimeSavedMs = 0;
end
if nTrial == 1 || ~isfield(input, 'didASave')
  input.didASave = false;
end

if input.subjectNum < 1
  disp('Must set subjectNum to save data - not saving');
  return
end

if ~exist(cs.dataPath, 'file')
  mkdir(cs.dataPath)
end
input.savedDataName = sprintf('%s/data-i%03d-%s.mat', ...
                              cs.dataPath, ...
                              input.subjectNum, datestr(now, 'yymmdd'));

%% compute elapsed time
if isempty(thisTimeMs)  % missing hold StartsMs; not sure when this happens
  doSave = false;
else
  savedDiffS = (thisTimeMs - input.lastTimeSavedMs)/1000;
  doSave = (mod(input.trialSinceReset, nTrialsToSkip) == 0 ...
            || savedDiffS > maxSToSkip); 
end



%% save vars to disk
if doSave
  % check and backup
  backup = {};

  % do we need a backup?
  if ~exist(input.savedDataName, 'file')
    backup = {};
  else
    ds = load(input.savedDataName);
  
    if ds.input.trialSinceReset < input.trialSinceReset  
      % fewer trials saved, save over it
      if isfield(ds, 'backup')
        backup = ds.backup;
      end
    else
      % move input to next backup field
      if isfield(ds, 'backup')
        backup = ds.backup;
      else
        backup = {};
      end
      backup{end+1} = ds.input;
    end
  end
  

  % save to disk
  if isempty(backup)
    save(input.savedDataName, 'input');
  else
    save(input.savedDataName, 'input', 'backup'); 
  end

  input.lastTimeSavedMs = thisTimeMs;
  input.didASave = true;
  disp(sprintf('Saved matlab variables to disk'));
end



