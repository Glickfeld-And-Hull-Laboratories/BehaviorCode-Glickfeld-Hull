function [x] = plotGlickfeldPerformance(subjMat)
% Rewritten Starting Dec 2014 to automatically detect experiment type and load all dates in
% the Data folder for a given mouse.
global dataStru
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = '~/Documents/MWorks/Data';
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
dataStru=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
    subjIx = subjNumList == j;
    desNames = fileNames(subjIx);
    nNames = length(desNames);
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        dataStru(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        
        %% hold time calculations
%        rawHoldTimes = cell2mat(ds.holdTimesMs);
%        medianHold = nanmedian(rawHoldTimes);
%        meanHold = nanmean(rawHoldTimes);
%        holdSTD = nanstd(double(rawHoldTimes));
%        uSTD = meanHold+holdSTD;
%%        tLSTD = meanHold-holdSTD;
%        if (meanHold-holdSTD)>0,
%          lSTD = tLSTD;
%        else
%          lSTD = 0;
%        end
        nTrials= length(ds.trialOutcomeCell);
        disp(ds.savedDataName)
        
        %% replace this part with hadc8/flashingStim2/lego sorting code
        if isfield(ds, 'tTotalReqHoldTimeMs'),
            dataStru(j).values.isLego(1, iN) = 0;
            %% do earlies analysis here
            earlyIx = strcmp(ds.trialOutcomeCell, 'failure');
            nEarly = sum(earlyIx);
            dataStru(j).values.nEarlies(1,iN)= nEarly;
            dataStru(j).values.perEarlies(1,iN)= nEarly/nTrials;
            dataStru(j).values.react(1,iN)= ds.reactTimeMs;
            
            if isfield(ds, 'targetOnTimeMs')
                dataStru(j).values.isLego(1, iN) = 0;
                dataStru(j).values.isHADC8(1, iN) = 0;
                dataStru(j).values.isFlashing(1, iN) = 1;
            else
                dataStru(j).values.isLego(1, iN) = 0;
                dataStru(j).values.isHADC8(1, iN) = 1;
                dataStru(j).values.isFlashing(1, iN) = 0;
            end
        elseif isfield(ds, 'delayTimeMs'),
            dataStru(j).values.isLego(1, iN) = 1;
            dataStru(j).values.isHADC8(1, iN) = 0;
            dataStru(j).values.isFlashing(1, iN) = 0;
            
            dataStru(j).values.leftIx{1, iN} = ds.tLeftTrial;
            
%            dataStru(j).values.minHold(1,iN) = ds.delayTimeMs+ds.tooFastTimeMs;
%            dataStru(j).values.maxHold(1,iN) = ds.delayTimeMs+ds.reactionTimeMs;
            
            %cancel out earlies analysis
            dataStru(j).values.nEarlies(1,iN)= NaN;
            dataStru(j).values.perEarlies(1,iN)= NaN;
            
            %% do separate incorrect analysis and avoid earlies
            
        end
        
        %% trial outcome calculations
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        nCorrect = sum(corrIx);
        igIx = strcmp(ds.trialOutcomeCell, 'ignore');
        nLate = sum(igIx);
       % corrHolds = rawHoldTimes(corrIx);
        
        perCorrect= nCorrect/nTrials;
        perLate= nLate/nTrials;
        %perEarly= nEarly/nTrials;
        incIx = strcmp(ds.trialOutcomeCell, 'incorrect')
        perInc= sum(incIx)/nTrials;
        dataStru(j).values.perInc(1,iN)= perInc;
        
%        dataStru(j).values.correctHoldSD(1, iN) = nanstd(double(corrHolds));
%        dataStru(j).values.correctHoldMean(1, iN) = nanmean(corrHolds);
        
        dataStru(j).values.nCorrects(1,iN)= nCorrect;
        dataStru(j).values.nLates(1,iN)= nLate;
        %dataStru(j).values.medianHold(1,iN)= medianHold;
        %dataStru(j).values.meanHold(1,iN)= meanHold;
        %dataStru(j).values.stdHold(1,iN)= holdSTD;
        %dataStru(j).values.uSTD(1,iN)= uSTD;
        %dataStru(j).values.lSTD(1,iN)= lSTD;
        dataStru(j).values.perCorrects(1,iN)= perCorrect;
        %dataStru(j).values.perEarlies(1,iN)= perEarly;
        dataStru(j).values.perLates(1,iN)= perLate;
        dataStru(j).values.ITI(1, iN)= ds.itiTimeMs;
        dataStru(j).values.nTrials(1, iN)= nTrials;
%        dataStru(j).values.holdTimesCell{iN} = histc(double(cell2mat(ds.holdTimesMs)), 0:100:5000);
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')
 

%%
subjMat = [513 514 518 520];
for i=1:length(subjMat),
    subPlotSz = {3,1};
    sub = subjMat(i);
    v = dataStru(sub).values;
    pCorr = v.perCorrects;
    pEarly = v.perEarlies;
    pLate = v.perLates; 
    pInc = v.perInc;
%    medH = v.medianHold;
%    uSTD = v.uSTD;
%    lSTD = v.lSTD;
    nTrials = v.nTrials;
    nCorr = v.nCorrects;
    nDays = length(v.perCorrects);
    xDays = 1:nDays;

%%    
    axH = subplot(subPlotSz{:}, 1);
    hold on
   
    plot(pCorr, '-kv', 'LineWidth', 2);
    plot(pEarly, '-cv', 'LineWidth', 2);
    plot(pLate, '-mv', 'LineWidth', 2);
    plot(pInc, '-gv', 'LineWidth', 2);
    ylabel('Percent');
    xlabel('Training Day');
    tName = strcat('Daily Trial Outcomes -- i', num2str(sub), '-- Generated:', datestr(today, 'mmmm dd yyyy'));
    title(tName)
    xlim([1 nDays])
    ylim([0 1]);
    %legend('Correct', 'Early', 'Late', 'No Release', 'Location', 'Best') ;
    grid on; 
%%
    axH = subplot(subPlotSz{:}, 2);
    hold on
    grid on; 
    
    lIx = cell2mat(v.leftIx);
    LCorr = pCorr(lIx);
    LLate = pLate(lIx); 
    LInc = pInc(lIx);
    
    plot(LCorr, '-kv', 'LineWidth', 2);
    plot(LLate, '-mv', 'LineWidth', 2);
    plot(LInc, '-gv', 'LineWidth', 2);
    ylabel('Left Percent');
    xlabel('Training Day');
    
%%    
    rIx = ~lIx;
    RCorr = pCorr(rIx);
    RLate = pLate(rIx); 
    RInc = pInc(rIx);
    
    plot(RCorr, '-kv', 'LineWidth', 2);
    plot(RLate, '-mv', 'LineWidth', 2);
    plot(RInc, '-gv', 'LineWidth', 2);
    ylabel('Right Percent');
    xlabel('Training Day');
    
 
    fName = strcat('Performance -- i', num2str(sub), ' -- Generated:', datestr(today, 'dd mmmm yyyy'));
    set(gcf, 'Name', fName);
    
    
    sName = strcat('~/Documents/MWorks/DailyPlots/plotPerformance-i', num2str(sub),'-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
         close
end