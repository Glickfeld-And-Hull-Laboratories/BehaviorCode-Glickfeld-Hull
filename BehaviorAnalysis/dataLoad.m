function [dataStru] = dataLoad(subjMat, dateMat)
%% Script to load data acquired on MWDataCentral into a structure (dataStruct) for analysis.
% To work, subjMat and datesMat must be 1xN matrices with YYMMDD format. Ex.[61 62 63 108 110
% 111...] or [130523 130524 130525...]

addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = '~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis/testdata';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
nSubjs= length(subjMat);

% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(4:end);

%% Data File Indexing and Pre-Processing
outR = regexp(fileNames, 'data-i([0-9]*)-([0-9]*).mat', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
dateStrList = str2double(outRM(:,2));
dataStru=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j= subjMat;
    subjIx = subjNumList == j;
    dateIx = ismember(dateStrList, dateMat);
    subjDateIx = subjIx & dateIx; 
    dataStru(j).dates= dateMat;
    desNames = fileNames(subjDateIx);
    nNames = length(desNames);
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        dataStru(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');    
        rawHoldTimes = cell2mat(ds.holdTimesMs);
        medianHold = nanmedian(rawHoldTimes);
        meanHold = nanmean(rawHoldTimes);
        nTrials= length(ds.trialOutcomeCell);
        
        nCorrects = sum(strcmp(ds.trialOutcomeCell, 'success'));
        nEarly = sum(strcmp(ds.trialOutcomeCell, 'early'));
        nLate = sum(strcmp(ds.trialOutcomeCell, 'late'));
        nNoRelease = sum(strcmp(ds.trialOutcomeCell, 'norelease'));       

        perCorrects= nCorrects/nTrials;
        perEarlies= nEarly/nTrials;
        perMisses= nLate/nTrials;
        perNoRelease = nNoRelease/nTrials;
  
        dataStru(j).values.nCorrects(1,iN)= nCorrects;
        dataStru(j).values.nEarlies(1,iN)= nEarly;
        dataStru(j).values.nMisses(1,iN)= nLate;
        dataStru(j).values.medianHold(1,iN)= medianHold;
        dataStru(j).values.meanHold(1,iN)= meanHold;
        dataStru(j).values.perCorrects(1,iN)= perCorrects;
        dataStru(j).values.perEarlies(1,iN)= perEarlies;
        dataStru(j).values.perMisses(1,iN)= perMisses;
        dataStru(j).values.delayTime(1, iN)= ds.delayTimeMs;
        dataStru(j).values.ITI(1, iN)= ds.itiTimeMs;
        dataStru(j).values.nTrials(1, iN)= nTrials;
       
       disp('Daily Data Loaded')  
    end
    disp('Subject Data Loaded')
end

beep
disp('Loading complete!')
