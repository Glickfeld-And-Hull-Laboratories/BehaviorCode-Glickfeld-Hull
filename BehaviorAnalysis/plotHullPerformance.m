function [x] = plotHullPerformance(subjMat, dateMat)
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
fileNames = fileNames(4:end);
%% Data File Indexing and Pre-Processing
outR = regexp(fileNames, 'data-i([0-9]*)-([0-9]*)-([0-9]*).mat', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
dateStrList = str2double(outRM(:,2));
hourStrList = str2double(outRM(:,3));
dataStru=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
    subjIx = subjNumList == j;
    if nargin<2,
        dateIx = dateStrList(subjIx);
    else
        dateIx = ismember(dateStrList, dateMat);
    end
    subjDateIx = subjIx & dateIx;
    %dataStru(j).dates= dateMat;
    desNames = fileNames(subjDateIx);
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
        
        nCorrect = sum(strcmp(ds.trialOutcomeCell, 'success'));
        nEarly = sum(strcmp(ds.trialOutcomeCell, 'failure'));
        nLate = sum(strcmp(ds.trialOutcomeCell, 'ignore'));
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
        dataStru(j).values.ITI(1, iN)= ds.itiTimeMs;
        dataStru(j).values.nTrials(1, iN)= nTrials;
        dataStru(j).values.holdTimesCell{iN} = histc(double(cell2mat(ds.holdTimesMs)), 0:100:5000);
        
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
        else
            % ASK EUNYOUNG HOW SHE DOES TIMING CONDITIONS W/O DNSC
        end
        
        
        disp('Daily Data Loaded')  
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')
       
for i=1:length(subjMat),
    subPlotSz = {3,1};
    sub = subjMat(i); 
    pCorr = dataStru(sub).values.perCorrects;
    pEarly = dataStru(sub).values.perEarlies;
    pLate = dataStru(sub).values.perLates;
    pNR = dataStru(sub).values.perNR;   
    medH = dataStru(sub).values.medianHold;
    uSTD = dataStru(sub).values.uSTD;
    lSTD = dataStru(sub).values.lSTD;
    nTrials = dataStru(sub).values.nTrials;
    nCorr = dataStru(sub).values.nCorrects;
    nDays = length(medH);
    xDays = 1:nDays;

    axH = subplot(subPlotSz{:}, 1);
    hold on
    plot(pCorr, '-kv', 'LineWidth', 2);
    plot(pEarly, '-cv', 'LineWidth', 2);
    plot(pLate, '-mv', 'LineWidth', 2);
    ylabel('Percent');
    xlabel('Training Day');
    tName = strcat('Trial Outcome Performance -- i', num2str(sub), '-- Generated:', datestr(today, 'dd mmmm yyyy'));
    title(tName);
    axis tight
    ylim([0 1]);
    %legend('Correct', 'Early', 'Late', 'No Release', 'Location', 'Best') ;
    grid on;

    axH = subplot(subPlotSz{:}, 3);
    hold on
    medPlot = plot(medH, 'k', 'LineWidth', 2);
%    winPlot = plot(winS, 'r', 'LineWidth', 2);
%    hTSPlot = plot(hTS, 'g', 'LineWidth', 2);
%    SD = jbfill(xDays, uSTD, lSTD, [0.25 0.25 1], [0.25 0.25 1], 1, 0.5);
    ylabel('Hold Time (ms)');
    xlabel('Training Day');
    axis tight
    grid on
    ylim([0 2500]);
    tName = strcat('Lever Hold Performance -- i', num2str(sub));%' -- Generated:', datestr(today, 'dd mmmm yyyy'));
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