addpath('C:\Users\andrew\Documents\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis');
% Data Path Formatting
global dataStruct
dataStruct = struct;
holds = {};
trs = {};

disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\andrew\Desktop\Analysis150520';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(3:end);

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
desNames = fileNames;
nNames = length(desNames);
for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        disp(tName) 
        dataStruct.values.dataName{iN} = tName;
        ds=mwLoadData(tName, 'max');
        
        hts = double(cell2mat(ds.holdTimesMs));
        hts(hts<200) = NaN;
        req = double(cell2mat(ds.tTotalReqHoldTimeMs));
        normalizedHolds = hts./req;
        trialNumbers = 1:length(hts);
        dataStruct.values.normalizedHolds{iN} = normalizedHolds;
        
        htStarts = ds.holdStartsMs;
        holdStarts = (cell2mat(htStarts) - htStarts{1})/60000;
        
        %% check if the first trial had randReq above 0 and doNoStimulusChange is not 1
        fTr = ds.firstTrConsts;
        if or(fTr.randReqHoldMaxMs<10, fTr.doNoStimulusChange==1),
            error(strcat('ERROR:', tName, ' does not begin with the Random Experiment. Rand: ', int2str(fTr.randReqHoldMaxMs), ', DoNoStimulusChange: ', int2str(fTr.randReqHoldMaxMs)))
        end
        
        %% get first trial constants
        dataStruct.values.startingRand{iN} = fTr.randReqHoldMaxMs;
        dataStruct.values.startingDNST{iN} = fTr.doNoStimulusChange;
        dataStruct.values.fixedStart{iN} = NaN;
        dataStruct.values.timingStart{iN} = NaN;
        
        %% look for variable changes
        changed = ds.constChangedOnTrial;
        changedOnTrialNumber = find(~cellfun(@isempty, changed));
        changes = changed(changedOnTrialNumber);
        for i=1:length(changes),
            j = length(changes{i});
            a = size(changes);
            for j=1:a(1),
                isRandChange = strcmp(changes{i}(j,1), 'randReqHoldMaxMs');
                if isRandChange==1,
                    dataStruct.values.randChanges{iN}(i) = changes{i}(j,2);
                    dataStruct.values.randChangeTr{iN}(i) = changedOnTrialNumber(i);
                    if cell2mat(changes{i}(j,2))<10,
                        dataStruct.values.fixedStart{iN} = changedOnTrialNumber(i);
                        %if changedOnTrialNumber(i)<50,
                        %   sprintf(strcat('Data file moved to fixed time in under 50 trials!'))
                        %end
                    end
                end
                isStimChange = strcmp(changes{i}(j,1), 'doNoStimulusChange');
                if isStimChange,
                    dataStruct.values.DNSCchanges{iN}(i) = changes{i}(j,2);
                    dataStruct.values.DNSCchangeTr{iN}(i) = changedOnTrialNumber(i);
                    if cell2mat(changes{i}(j,1))==1,
                        dataStruct.values.timingStart{iN} = changedOnTrialNumber(i);
                    end
                end
            end
        end
        
        %% crop holds and trials based on variable changes
        if isnan(dataStruct.values.fixedStart{iN}),
             dataStruct.values.randomHolds{iN} = normalizedHolds;
             dataStruct.values.randomTrs{iN} = trialNumbers;
             dataStruct.values.randomTime{iN} = holdStarts;
        else
            dataStruct.values.randomHolds{iN} = normalizedHolds(1:(double(dataStruct.values.fixedStart{iN})-1));
            dataStruct.values.randomTrs{iN} = trialNumbers(1:(double(dataStruct.values.fixedStart{iN})-1));
            dataStruct.values.randomTime{iN} = holdStarts(1:(double(dataStruct.values.fixedStart{iN})-1));
        end
        if isnan(dataStruct.values.timingStart{iN}),
            dataStruct.values.fixedHolds{iN} = normalizedHolds(dataStruct.values.fixedStart{iN}:length(normalizedHolds));
            dataStruct.values.fixedTrs{iN} = trialNumbers(dataStruct.values.fixedStart{iN}:length(normalizedHolds));
            dataStruct.values.fixedTime{iN} = holdStarts(dataStruct.values.fixedStart{iN}:length(normalizedHolds));
        else
            dataStruct.values.fixedHolds{iN} = normalizedHolds(dataStruct.values.fixedStart{iN}:dataStruct.values.timingStart{iN});
            dataStruct.values.fixedTrs{iN} = trialNumbers(dataStruct.values.fixedStart{iN}:dataStruct.values.timingStart{iN});
            dataStruct.values.fixedTime{iN} = holdStarts(dataStruct.values.fixedStart{iN}:dataStruct.values.timingStart{iN});
        end
       disp('Daily Data Loaded')
end
%%
[sortedRandomTrs ix] = sort(cell2mat_padded(dataStruct.values.randomTrs));
randomHolds = cell2mat_padded(dataStruct.values.randomHolds);
sortedRandomHolds = randomHolds(ix);


holdMeans = {};
holdSTDs = {};
x = 0:10:1990;
y = 10:10:2000;
i=1;
for i = 1:length(x),
    
    iX = x(i)<sortedRandomTrs & sortedRandomTrs<y(i);
    trIx = sortedRandomHolds(iX);
    holdSTDs{i} = nanstd(double(trIx));
    holdMeans{i} = nanmean(double(trIx));
end
figure;
hold on
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
plot([0 1200], [1 1], 'k')
grid on
xlim([0 200])
ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)')
%ylabel('Mean Reaction Time (ms)')
ylim([0 3])
title('Random-to-Fixed Time Experiment Summary: Random Trials')
xlabel('Trials')

sName = 'C:\Users\andrew\Desktop\randomSummary-trials.pdf';
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
    close

%%
[sortedFixedTrs ix] = sort(cell2mat_padded(dataStruct.values.fixedTrs));
fixedHolds = cell2mat_padded(dataStruct.values.fixedHolds);
sortedFixedHolds = fixedHolds(ix);


holdMeans = {};
holdSTDs = {};
x = 0:10:1990;
y = 10:10:2000;
i=1;
for i = 1:length(x),
    
    iX = x(i)<sortedFixedTrs & sortedFixedTrs<y(i);
    trIx = sortedFixedHolds(iX);
    holdSTDs{i} = nanstd(double(trIx));
    holdMeans{i} = nanmean(double(trIx));
end
figure;
hold on
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
plot([0 1200], [1 1], 'k')
grid on
xlim([0 1000])
ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)')
%ylabel('Mean Reaction Time (ms)')
ylim([0 3])
title('Fixed-to-Random Experiment Summary: Fixed Time')
xlabel('Trials')
sName = 'C:\Users\andrew\Desktop\fixedSummary-trials.pdf';
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
    close

%%
[sortedRandomTime ix] = sort(cell2mat_padded(dataStruct.values.randomTime));
randomHolds = cell2mat_padded(dataStruct.values.randomHolds);
sortedRandomHolds = randomHolds(ix);

holdMeans = {};
holdSTDs = {};
x = 0:2:118;
y = 2:2:120;
for i = 1:length(x),
    
    iX = x(i)<sortedRandomTime & sortedRandomTime<y(i);
    trs = sortedRandomHolds(iX);
    holdSTDs{i} = nanstd(double(trs));
    holdMeans{i} = nanmean(double(trs));
end
figure;
hold on
plot([0 1200], [1 1], 'k');
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
grid on
ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)');
%ylabel('Mean Reaction Time (ms)')
ylim([0 3])
xlim([1 30])
title('Fixed-to-Random Experiment Summary: Random Time')
xlabel('Minutes in Training')
sName = 'C:\Users\andrew\Desktop\randomSummary-time.pdf';
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
    close

%%
[sortedFixedTime ix] = sort(cell2mat_padded(dataStruct.values.fixedTime));
fixedHolds = cell2mat_padded(dataStruct.values.fixedHolds);
sortedFixedHolds = fixedHolds(ix);

holdMeans = {};
holdSTDs = {};
x = 0:2:118;
y = 2:2:120;
for i = 1:length(x),
    
    iX = x(i)<sortedFixedTime & sortedFixedTime<y(i);
    trs = sortedFixedHolds(iX);
    holdSTDs{i} = nanstd(double(trs));
    holdMeans{i} = nanmean(double(trs));
end
figure;
hold on
plot([0 1200], [1 1], 'k');
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
grid on
ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)');
%ylabel('Mean Reaction Time (ms)')
ylim([0 3])
xlim([1 100])
title('Fixed-to-Random Experiment Summary: Fixed Time')
xlabel('Minutes in Training')
sName = 'C:\Users\andrew\Desktop\fixedSummary-time.pdf';
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
    close

