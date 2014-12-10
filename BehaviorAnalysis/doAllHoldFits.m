function [dataStru outR] = doAllHoldFits(subjMat, dateMat)
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
colorCell = {[0 0 1], ... % blue
             [0 1 0], ... & green
             [0.8 0.8 0], ... % yellow
             [1 0.5 0], ... % orange
             [1 0 0], ... % red
             [0.8 0 0.8]}; % purple
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(3:end);
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
    nNames
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        dataStru(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        [a b] = fitDistToHolds(ds)
        if iN>length(colorCell);
            t = iN - length(colorCell);
        else
            t = iN;
        end
        set(a, 'Color', colorCell{t});
        set(b, 'Color', colorCell{t});
        hold on
        
       disp('Daily Data Loaded')  
    end
    disp('Subject Data Loaded')
end

beep
disp('Loading complete!')

%% Plotting Functions Here
%plotPerformance(subjMat)