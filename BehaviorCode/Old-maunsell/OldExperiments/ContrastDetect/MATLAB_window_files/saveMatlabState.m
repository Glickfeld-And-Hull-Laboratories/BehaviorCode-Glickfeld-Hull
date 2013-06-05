function input = saveMatlabState(data_struct, input)

% compute whether to do export on this trial
nTrialsToSkip = 10;
maxSToSkip = 120;

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
  error('Must set subjectNum to save data');
end

input.savedDataName = sprintf('~/Library/Application Support/MWorks/DataFiles/data-i%03d-%s.mat', ...
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
  if exist(input.savedDataName)
    ds = load(input.savedDataName);
    if isfield(ds, 'backup')
      backup = ds.backup;
    end
    if ~input.didASave  % append old data only if we haven't saved since reset
      backup{end+1} = ds;
    end  
  end
  
  if isempty(backup)
    save(input.savedDataName, 'input');
  else
    save(input.savedDataName, 'input', 'backup'); 
  end

  input.lastTimeSavedMs = thisTimeMs;
  input.didASave = true;
  disp(sprintf('Saved matlab variables to disk'));
end



