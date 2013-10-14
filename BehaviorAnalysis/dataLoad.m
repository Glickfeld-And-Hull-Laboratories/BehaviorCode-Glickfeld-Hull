function [dataStru] = dataLoad(subjMat, dateMat)
%% Script to load data acquired on MWDataCentral into a structure (dataStruct) for analysis.
% To work, subjMat and datesMat must be 1xN matrices with YYMMDD format. Ex.[61 62 63 108 110
% 111...] or [130523 130524 130525...]

global dataStru
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = '~/Documents/MWorks/Data';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
nSubjs= length(subjMat);
if nargin<2,
    a = today;
    d = (735472:[today]);
    dateMat = datestr(d, 'yymmdd');
    dateMat = str2num(dateMat);
end

% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(11:end);
%% Data File Indexing and Pre-Processing
outR = regexp(fileNames, 'data-i([0-9]*)-([0-9]*).mat', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
dateStrList = str2double(outRM(:,2));
dataStru=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
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
        
        % This bit makes windowStart calculating code compatible with early window width code.
        if isfield(ds, 'preRewardWindowMs'),
            dataStru(j).values.windowStart(1, iN) = ds.delayTimeMs - ds.preRewardWindowMs;
        elseif isfield(ds, 'rewardWindowWidthMs'),
            rewWindowStart = ds.delayTimeMs - (ds.rewardWindowWidthMs./2);
            if rewWindowStart > 0;
                dataStru(j).values.windowStart(1, iN) = rewWindowStart;
            else
                dataStru(j).values.windowStart(1, iN) = 0;
            end
        end
        
        rawHoldTimes = cell2mat(ds.holdTimesMs);
        medianHold = nanmedian(rawHoldTimes);
        meanHold = nanmean(rawHoldTimes);
        holdSTD = nanstd(double(rawHoldTimes));
        uSTD = meanHold+holdSTD;
        tLSTD = meanHold-holdSTD;
        if (meanHold-holdSTD)>0,
          lSTD = tLSTD;
        else
          lSTD = 0;
        end
        nTrials= length(ds.trialOutcomeCell);
        
        nCorrect = sum(strcmp(ds.trialOutcomeCell, 'success'));
        nEarly = sum(strcmp(ds.trialOutcomeCell, 'early'));
        nLate = sum(strcmp(ds.trialOutcomeCell, 'late'));
        nNR = sum(strcmp(ds.trialOutcomeCell, 'norelease'));       

        perCorrect= nCorrect/nTrials;
        perEarly= nEarly/nTrials;
        perLate= nLate/nTrials;
        perNR = nNR/nTrials;
  
        dataStru(j).values.nCorrects(1,iN)= nCorrect;
        dataStru(j).values.nEarlies(1,iN)= nEarly;
        dataStru(j).values.nLates(1,iN)= nLate;
        dataStru(j).values.nNR(1,iN)= nNR;
        dataStru(j).values.medianHold(1,iN)= medianHold;
        dataStru(j).values.meanHold(1,iN)= meanHold;
        dataStru(j).values.stdHold(1,iN)= holdSTD;
        dataStru(j).values.uSTD(1,iN)= uSTD;
        dataStru(j).values.lSTD(1,iN)= lSTD;
        dataStru(j).values.perCorrects(1,iN)= perCorrect;
        dataStru(j).values.perEarlies(1,iN)= perEarly;
        dataStru(j).values.perLates(1,iN)= perLate;
        dataStru(j).values.perNR(1,iN)= perNR;
        dataStru(j).values.delayTime(1, iN)= ds.delayTimeMs;
        dataStru(j).values.ITI(1, iN)= ds.itiTimeMs;
        dataStru(j).values.nTrials(1, iN)= nTrials;
       
       disp('Daily Data Loaded')  
    end
    disp('Subject Data Loaded')
end

beep
disp('Loading complete!')

%% Plotting Functions Here
plotPerformance(subjMat)

