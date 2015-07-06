function [x] = plotGlickfeldPerformance(subjMat)
% Rewritten Starting Dec 2014 to automatically detect experiment type and load all dates in
% the Data folder for a given mouse.
global dataStruct
if ispc,
    disp('I AM A PC')
    addpath('C:\Users\jake\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis\');
elseif isunix,
    disp('I AM A MAC OR LINUX COMPUTER')
    addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
end
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\jake\TempData\BxDatabase';
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
    nNames = length(desNames)
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        dataStruct(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        
        nTrials = length(ds.trialOutcomeCell);
        dataStruct(j).values.nTrials{1,iN} = nTrials;
        %% trial outcome calculations
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        dataStruct(j).values.perCorrect{1,iN} = sum(corrIx)/nTrials;
        igIx = strcmp(ds.trialOutcomeCell, 'ignore');
        dataStruct(j).values.perIgnore{1,iN} = sum(igIx)/nTrials;
        incIx = strcmp(ds.trialOutcomeCell, 'incorrect');
        dataStruct(j).values.perIncorrect{1,iN} = sum(incIx)/nTrials;
        earlyIx = strcmp(ds.trialOutcomeCell, 'failure');
        dataStruct(j).values.perEarly{1,iN} = sum(earlyIx)/nTrials;
        a = mean(cell2mat(ds.tTotalReqHoldTimeMs));
        b = mean(cell2mat(ds.holdTimesMs));
        dataStruct(j).values.meanReqHold{1,iN} = a;
        dataStruct(j).values.meanActualHold{1,iN} = b;
        reacts = cell2mat(ds.reactTimesMs);
        dataStruct(j).values.meanReact{1,iN} = mean(reacts(corrIx));
        
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')
 

%%
for i=subjMat,
    subPlotSz = {3,1};
    v = dataStruct(i).values;
    pCorr = cell2mat(v.perCorrect);
 %   pEarly = v.perEarlies;
    pLate = cell2mat(v.perIgnore); 
    pEarly = cell2mat(v.perEarly);
%    medH = v.medianHold;
%    uSTD = v.uSTD;
%    lSTD = v.lSTD;
    nDays = length(v.perCorrect);
    xDays = 1:nDays;

%%  

    axH = subplot(subPlotSz{:}, 1);
    hold on

    plot(pCorr, '-kv', 'LineWidth', 2);
    plot(pLate, '-mv', 'LineWidth', 2);
    plot(pEarly, '-cv', 'LineWidth', 2);
    ylabel('Percent');
    xlabel('Training Day');
    tName = strcat('Daily Trial Outcomes -- i', num2str(i), '-- Generated:', datestr(now, 'mmmm dd yyyy'));
    title(tName)
    xlim([1 nDays])
    ylim([0 1]);
    grid on
    
    axH = subplot(subPlotSz{:}, 2);
    hold on
    
    plot(cell2mat(v.meanReqHold), '-kv', 'LineWidth', 2);
    plot(cell2mat(v.meanActualHold), '-rv', 'LineWidth', 2);
    ylabel('Hold Time (ms)');
    xlabel('Training Day');
    tName = strcat('Required (black) vs Actual Holds (red)');
    title(tName)
    xlim([1 nDays])
    %ylim([0 1]);
    grid on

    axH = subplot(subPlotSz{:}, 3);
    meanReact = v.meanReact;
    nTrials = v.nTrials;
    [axH, a, b] = plotyy(1:length(nTrials), cell2mat(nTrials), 1:length(meanReact), cell2mat(meanReact));
    xlim([1 length(nTrials)])
    ylabel(axH(1), 'Number of Trials');
    ylabel(axH(2), 'Mean Reaction Time (ms)');
    title('Number of Trials and Mean Reaction Times');
    grid on
    
    fName = strcat('Performance -- i', num2str(i), ' -- Generated:', datestr(now, 'dd mmmm yyyy'));
    set(gcf, 'Name', fName);
    
    
    sName = strcat('C:\Users\jake\TempData\', num2str(i),'-', datestr(now, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
    close
end
end