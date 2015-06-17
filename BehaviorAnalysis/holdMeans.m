function [dataStruct outR] = holdMeans(subjMat)
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\andrew\Desktop\behavData';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
nSubjs= length(subjMat);
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(3:end);
%% Data File Indexing and Pre-Processing
outR = regexp(fileNames, 'data-i([0-9]*)-([0-9]*)([???]*)', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
dateStrList = str2double(outRM(:,2));
dataStru=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
    subjIx = subjNumList == j;
    subjDateIx = subjIx;
    %dataStru(j).dates= dateMat;
    desNames = fileNames(subjDateIx);
    nNames = length(desNames);
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN})
        dataStruct(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        rt = cell2mat(ds.reactTimesMs);
        dataStruct(j).values.reactMeans{iN} = mean(rt);
        dataStruct(j).values.reactSD{iN} = std(double(rt));
        if ds.randReqHoldMaxMs>10,
            dataStruct(j).values.fixedDay(iN) = 0;
        else
            dataStruct(j).values.fixedDay(iN) = 1;
        end
        to = ds.trialOutcomeCell;
        corrIx = strcmp(to, 'success');
        dataStruct(j).values.corrReactMeans{iN} = mean(rt(corrIx));
        dataStruct(j).values.corrReactSD{iN} = std(double(rt(corrIx)));
        
       disp('Daily Data Loaded')  
    end
    disp('Subject Data Loaded')
    a = [1:length(dataStruct(j).values.reactMeans)]
    [b c d] = plotyy(a, cell2mat(dataStruct(j).values.reactSD), a, dataStruct(j).values.fixedDay)
    set(c,'Marker','o')
    title('Standard Deviation Data')
    xlabel('Training Day')
    ylabel('Reaction Time Standard Deviation')
    
    figure;
    [b c d] = plotyy(a, cell2mat(dataStruct(j).values.corrReactSD), a, dataStruct(j).values.fixedDay)
    set(c,'Marker','o')
    title('Standard Deviation Data: Corrects Only')
    xlabel('Training Day')
    ylabel('Reaction Time Standard Deviation')
    
    
    
end

beep
disp('Loading complete!')
