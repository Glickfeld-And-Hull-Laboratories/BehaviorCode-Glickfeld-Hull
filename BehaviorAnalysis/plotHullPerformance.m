function [x] = plotHullPerformance(subjMat)
% Rewritten December 16th 2014 to automatically detect experiment type and load all dates in
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
        disp(ds.savedDataName)
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        nCorrect = sum(corrIx);
        corrHolds = rawHoldTimes(corrIx);
        
        dataStru(j).values.correctHoldSD(1, iN) = nanstd(double(corrHolds));
        dataStru(j).values.correctHoldMean(1, iN) = nanmean(corrHolds);
        
        nLate = sum(strcmp(ds.trialOutcomeCell, 'ignore'));
        nNR = sum(strcmp(ds.trialOutcomeCell, 'norelease'));
        perCorrect= nCorrect/nTrials;
        perNR = nNR/nTrials;
        
        if isfield(ds, 'tTotalReqHoldTimeMs'),
            dataStru(j).values.minReq(1, iN) = min(cell2mat(ds.tTotalReqHoldTimeMs))+ds.tooFastTimeMs;
            dataStru(j).values.maxReq(1, iN) = max(cell2mat(ds.tTotalReqHoldTimeMs))+ds.reactTimeMs;
            nLate = sum(strcmp(ds.trialOutcomeCell, 'ignore'));
            nEarly = sum(strcmp(ds.trialOutcomeCell, 'failure'));
        elseif isfield(ds, 'delayTimeMs'),
            dataStru(j).values.minReq(1, iN) = min(cell2mat(ds.tDelayTimeMs)-ds.preRewardWindowMs);
            dataStru(j).values.maxReq(1, iN) = max(cell2mat(ds.tDelayTimeMs)+ds.postRewardWindowMs);
            nLate = sum(strcmp(ds.trialOutcomeCell, 'late'));
            nEarly = sum(strcmp(ds.trialOutcomeCell, 'early'));
        end
        perLate= nLate/nTrials;
        perEarly= nEarly/nTrials;
        
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
        dataStru(j).values.ITI(1, iN)= ds.itiTimeMs;
        dataStru(j).values.nTrials(1, iN)= nTrials;
        dataStru(j).values.holdTimesCell{iN} = histc(double(cell2mat(ds.holdTimesMs)), 0:100:5000);
        if isfield(ds, 'tTotalReqHoldTimeMs')
            if isfield(ds, 'doNoStimulusChange')
                if ds.doNoStimulusChange == 1,
                    dataStru(j).values.timing(1, iN)= 1;
                    dataStru(j).values.fixedCue(1, iN)= 0;
                    dataStru(j).values.randCue(1, iN)= 0;
                % Then calculate things for timing conditions
                elseif ds.doNoStimulusChange == 0,
                    dataStru(j).values.timing(1, iN)= 0;
                    if ds.randReqHoldMaxMs > 5,
                        % Calculate things for random cued conditions
                        dataStru(j).values.randCue(1, iN)= 1;
                        dataStru(j).values.fixedCue(1, iN)= 0;
                    else
                        % Calculate things for fixed cued conditions
                        dataStru(j).values.randCue(1, iN)= 0;
                        dataStru(j).values.fixedCue(1, iN)= 1;
                    end
                end
            else %disp('ERROR IN CATEGORIZATION: HADC8 without doNoStimulusChange')
                if ds.randReqHoldMaxMs > 5,
                    dataStru(j).values.timing(1, iN)= 0;
                    dataStru(j).values.fixedCue(1, iN)= 0;
                    dataStru(j).values.randCue(1, iN)= 1;
                else
                    dataStru(j).values.timing(1, iN)= 0;
                    dataStru(j).values.fixedCue(1, iN)= 1;
                    dataStru(j).values.randCue(1, iN)= 0;
                end
            end
        elseif isfield(ds, 'delayTimeMs'),
            dataStru(j).values.timing(1, iN)= 1;
            dataStru(j).values.randCue(1, iN)= 0;
            dataStru(j).values.fixedCue(1, iN)= 0;
        end
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')
 

%%
for i=1:length(subjMat),
    subPlotSz = {3,1};
    sub = subjMat(i);
    v = dataStru(sub).values;
    pCorr = v.perCorrects;
    pEarly = v.perEarlies;
    pLate = v.perLates;
    pNR = v.perNR;   
    medH = v.medianHold;
    uSTD = v.uSTD;
    lSTD = v.lSTD;
    nTrials = v.nTrials;
    nCorr = v.nCorrects;
    nDays = length(medH);
    xDays = 1:nDays;

%%    
    axH = subplot(subPlotSz{:}, 1);
    
    cueIx = v.fixedCue|v.randCue;
    hold on
    cueCorr = pCorr;
    cueCorr(~cueIx) = NaN;
    cueEarly = pEarly;
    cueEarly(~cueIx) = NaN;
    cueLate = pLate;
    cueLate(~cueIx) = NaN;
    if sum(v.perNR>0),
        timingNR = v.perNR;
        timingNR(~cueIx) = NaN;
        plot(timingNR, '-v', 'Color', [1 0.6 0])
    end
   
    plot(cueCorr, '-kv', 'LineWidth', 2);
    plot(cueEarly, '-cv', 'LineWidth', 2);
    plot(cueLate, '-mv', 'LineWidth', 2);
    ylabel('Percent (Cued)');
    xlabel('Training Day');
    tName = strcat('Trial Outcomes for Cued Conditions -- i', num2str(sub), '-- Generated:', datestr(today, 'mmmm dd yyyy'));
    title(tName)
    xlim([1 nDays])
    ylim([0 1]);
    %legend('Correct', 'Early', 'Late', 'No Release', 'Location', 'Best') ;
    grid on; 
%%
    axH = subplot(subPlotSz{:}, 2);
    
    tIx = logical(v.timing)
    hold on
    timingCorr = pCorr;
    timingCorr(~tIx) = NaN;
    timingEarly = pEarly;
    timingEarly(~tIx) = NaN;
    timingLate = pLate;
    timingLate(~tIx) = NaN;
    if sum(v.perNR>0),
        timingNR = v.perNR;
        timingNR(~tIx) = NaN;
        plot(timingNR, '-v', 'Color', [1 0.6 0])
    end
        
   
    plot(timingCorr, '-kv', 'LineWidth', 2);
    plot(timingEarly, '-cv', 'LineWidth', 2);
    plot(timingLate, '-mv', 'LineWidth', 2);
    ylabel('Percent (Timing)');
    xlabel('Training Day');
    tName = 'Trial Outcomes for Timing Conditions';
    title(tName);
    xlim([1 nDays])
    ylim([0 1]);
    %legend('Correct', 'Early', 'Late', 'No Release', 'Location', 'Best') ;
    grid on;

%%
    axH = subplot(subPlotSz{:}, 3);
    hold on
    medPlot = plot(medH, 'k', 'LineWidth', 2);
    [ax] = jbfill(xDays, v.maxReq, v.minReq, [0 0 1], [0 0 1], 1, [0.5]);
    anystack(ax, 'bottom');
%    winPlot = plot(winS, 'r', 'LineWidth', 2);
%    hTSPlot = plot(hTS, 'g', 'LineWidth', 2);
    ylabel('Hold Time (ms)');
    xlabel('Training Day');
    axis tight
    grid on
    m = max(v.maxReq)
    if m>5000,
        ylim([0 5000]);
    else
        ylim([0 m]);
    end
    tName = 'Median Hold Times and Reward Window Sizes';
    title(tName);
    %legend('\pm 1 SD', 'Median Hold Time', 'Reward Window Start', 'Location', 'Best');
    
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