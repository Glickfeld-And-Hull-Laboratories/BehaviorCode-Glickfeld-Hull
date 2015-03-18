function [x] = plotGlickfeldPerformance(subjMat)
subjMat = [519]
% Rewritten Starting Dec 2014 to automatically detect experiment type and load all dates in
% the Data folder for a given mouse.
global dataStru
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\andrew\Desktop\BehavData';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
nSubjs= length(subjMat);
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(5:end);
%% Data File Indexing and Pre-Processing
outR = regexp(fileNames, 'data-i([0-9]*)-(??????????????', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
%dateStrList = str2double(outRM(:,2));
%hourStrList = str2double(outRM(:,3));
dataStruct=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
    subjIx = subjNumList == j;
    desNames = fileNames(subjIx);
    nNames = length(desNames);
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        dataStruct(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        
        nTrials = length(ds.trialOutcomeCell);
        dataStruct(j).values.nTrials = nTrials;
        %% trial outcome calculations
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        dataStruct(j).values.perCorrect{1,iN} = sum(corrIx)/nTrials;
        igIx = strcmp(ds.trialOutcomeCell, 'ignore');
        dataStruct(j).values.perIgnore{1,iN} = sum(igIx)/nTrials;
        incIx = strcmp(ds.trialOutcomeCell, 'incorrect');
        dataStruct(j).values.perIncorrect{1,iN} = sum(incIx)/nTrials;
        
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')
 

%%
for i=1:length(subjMat),
    subPlotSz = {3,1};
    sub = subjMat(i);
    v = dataStruct(519).values;
    pCorr = cell2mat(v.perCorrect);
 %   pEarly = v.perEarlies;
    pLate = cell2mat(v.perIgnore); 
    pInc = cell2mat(v.perIncorrect);
%    medH = v.medianHold;
%    uSTD = v.uSTD;
%    lSTD = v.lSTD;
    nDays = length(v.perCorrect);
    xDays = 1:nDays;

%%  

    axH = subplot(subPlotSz{:}, 1:3);
    hold on

    plot(pCorr, '-kv', 'LineWidth', 2);
   % plot(pEarly, '-cv', 'LineWidth', 2);
    plot(pLate, '-mv', 'LineWidth', 2);
    plot(pInc, '-gv', 'LineWidth', 2);
    ylabel('Percent');
    xlabel('Training Day');
    tName = strcat('Daily Trial Outcomes -- i', num2str(sub), '-- Generated:', datestr(today, 'mmmm dd yyyy'));
    title(tName)
    xlim([1 nDays])
    ylim([0 1]);
    %legend('Correct', 'Early', 'Late', 'No Release', 'Location', 'Best')

        disp('Fig 1 failed');
    
    fName = strcat('Performance -- i', num2str(sub), ' -- Generated:', datestr(today, 'dd mmmm yyyy'));
    set(gcf, 'Name', fName);
    
    
    sName = strcat('C:\Users\andrew\Desktop\BehavData\', num2str(sub),'-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
end
end