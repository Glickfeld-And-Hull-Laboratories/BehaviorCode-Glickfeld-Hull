clear all
close all

%% Update behavior indices sheet

updateBehaviorIndices = 0 %set to 1 if you want to update the excel spreadsheet before running, set to 0 otherwise


if updateBehaviorIndices == 1;
    eval('UpdateBehaviorIndices')
end

%% Specify mice

allmice = ["i2198"]
modelist = ["train"] 
%allmice = ["i2182", "i2179", "i2174", "i2186"]
%modelist = ["test", "test", "test", "train"]
saveImage = 1

for i = 1:length(allmice)

%% Specify variables

mouse = allmice(i)
mode = modelist(i) %test or train
interval = 5 %number of days worth of data you want to see
endChop = 0 %number of days to cut off the end (useful for debugging)
thresholdHitRate = 0; %set threshold for hit rate on easiest trials, sessions with hit rate below threshold won't be included
requireBlock2 = 0; %set to 1 if you only want to inlude sessions with block 2 enabled, otherwise set to 0 to include all sessions


%% Load behavioral data

%Import excel spreadsheet with file info
%(note that you may have to update the columns you use as you go using the "Range" parameter)
excel = readtable("home\LindseyW\MATLAB\Behavior_analysis\BehaviorIndices.xlsx", Range="A:D", Sheet = 1, ReadVariableNames = true);

miceList = excel{:,2}; %isolates the column of the excel table where the mouse name is stored
goodSessions = find(string(miceList) == mouse); %filters out just the rows with the mouse we want
filenames = string(excel{goodSessions,1}); %gives us the file names associated with those filtered out sessions
recentSessions = filenames(length(filenames)-(interval-1):length(filenames)-endChop); %loads just the most recent n sessions

%Manual file name entry (COMMENT OFF WHEN NOT IN USE)
%recentSessions = ["data-i2167-241211-1203.mat", "data-i2167-241213-1051.mat", "data-i2167-241214-1130.mat", "data-i2167-241215-1035.mat"]'


hitRateSummary = zeros(3, length(recentSessions));

%%testing mode
if mode == "test"
    dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
    for ifile = 1:length(recentSessions); %COME BACK AND CHANGE TO 1
       try %To catch errors caused by corrupt files
            load(fullfile(dir,recentSessions(ifile)));
       catch
            disp("One or more files were corrupt and did not load")
            continue
        end 
           
    
        %pull put performance variables
        cons = celleqel2mat_padded(input.tGratingContrast);
        conlist = unique(cons); %list of contrasts presented in the session
        maxCon = max(conlist); %finds the max contrast (easiest contrast to see)
        SIx = strcmp(input.trialOutcomeCell,'success'); %Index of successes (ones)
        MIx = strcmp(input.trialOutcomeCell,'ignore'); %Index of misses (ones)
        FAIx = strcmp(input.trialOutcomeCell,'failure'); %Index of false alarms (ones)
    
    
        %checks hit rate for trials with the highest contrast
        easyTrials = find(cons == maxCon);
        hitOnEasiest = length(intersect(find(SIx),easyTrials));
        missOnEasiest = length(intersect(find(MIx),easyTrials));
        hitRateOnEasiest = hitOnEasiest/(hitOnEasiest + missOnEasiest);
        hitRateSummary(1, ifile) = hitOnEasiest;
        hitRateSummary(2, ifile) = missOnEasiest;
        hitRateSummary(3, ifile) = hitRateOnEasiest;
        
    
        %add the sessions that meet your criteria to a temporary structure
        if requireBlock2 == 1;
            if (block2On == 1) && (hitRateOnEasiest >= thresholdHitRate);
                temp(ifile) = input; 
            else
                continue
            end
        elseif requireBlock2 == 0;
            if hitRateOnEasiest >= thresholdHitRate;
                try %Again, to get around corrupt files
                    temp(ifile) = input; %Structs from imaging sessions are different from rgular structs, and this won't join them
                catch
                    continue
                end 
                    %}
            else
                continue
            end
        end
    
    
    
    %{
    %Concatenate the data in the temporary 1 x n structure 
    % into a single 1 x 1 with all of your selected trials
    all_data = concatenateDataBlocks(temp);
    %}
    
    %% plot online performance
    
    
        % Online performance
        trialNum = linspace(1, length(input.tGratingContrast), length(input.tGratingContrast));
        nTrial = length(trialNum);
    
      
    
        dateText = string(input.saveTime);
        fixedTime = string(input.fixedReqHoldTimeMs);
        randTime = string(input.randReqHoldMaxMs);
        tooFastTime = string(input.tooFastTimeMs);
        parameterText = join(["Fixed =", fixedTime, "Random =", randTime, "TooFast =", tooFastTime]);
    
        figure(i);
        sgtitle(mouse);
        subplot(3,interval,ifile)
        hold on;
        if ifile <= 2
        perf1 = plot(smooth(SIx, ceil(nTrial/10), 'lowess'));
        ylim([0 1]);
       
        perf2 = plot(smooth(SIx, nTrial, 'lowess'));
        set(perf2, 'Color', 'r', ...
                'Linewidth', 3)
         ylim([0 1]);
        
        perf3 = plot(smooth(SIx, 100, 'lowess'));
        set(perf3, 'Color', 'k', ...
            'LineWidth', 2);
         ylim([0 1]);
         
        perf4 = plot(smooth(MIx, 100, 'lowess'));
        set(perf4, 'Color', 'm', ...
                'LineWidth', 2, ...
                'LineStyle', '-.');
         ylim([0 1]);
       
        perf5 = plot(smooth(FAIx, 100, 'lowess'));
        set(perf5, 'Color', [0 0.6 0], ...
                'LineWidth', 1, ...
                'LineStyle', '-.');
        ylim([0 1]);

        xlabel('Trial Number')
        ylabel('Rate')
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
        title(dateText)
        subtitle(parameterText, "FontSize", 8)
        else

        perf1 = plot(smooth(SIx, ceil(nTrial/10), 'lowess'));
        set(perf1, 'DisplayName', 'Hit 1');
         ylim([0 1]);
       
        perf2 = plot(smooth(SIx, nTrial, 'lowess'));
        set(perf2, 'Color', 'r', ...
                'Linewidth', 3, 'DisplayName', 'Hit 2')
         ylim([0 1]);
        
        perf3 = plot(smooth(SIx, 100, 'lowess'));
        set(perf3, 'Color', 'k', ...
            'LineWidth', 2, 'DisplayName', 'Hit 3');
         ylim([0 1]);
         
        perf4 = plot(smooth(MIx, 100, 'lowess'));
        set(perf4, 'Color', 'm', ...
                'LineWidth', 2, ...
                'LineStyle', '-.', 'DisplayName', 'Miss');
         ylim([0 1]);
       
        perf5 = plot(smooth(FAIx, 100, 'lowess'));
        set(perf5, 'Color', [0 0.6 0], ...
                'LineWidth', 1, ...
                'LineStyle', '-.', 'DisplayName', 'FA');
        ylim([0 1]);

        xlabel('Trial Number')
        ylabel('Rate')
        title(dateText)
        subtitle(parameterText, "FontSize", 8)
    
        Lgnd = legend('show');
        Lgnd.Position(1) = 0.02;
        Lgnd.Position(2) = 0.73;
            %legend('Hit', 'Hit', 'Hit(Mid)', 'Miss', 'FA', 'Location','eastoutside')
        end
    
    
    
    %% Plot psychometric curves
        
    
        % Extract trial outcomes
        tGratingContrast = celleqel2mat_padded(input.tGratingContrast); %contrast for each trial
        
        b2IxRaw = celleqel2mat_padded(input.tBlock2TrialNumber); %IDs trial as block1 (zeros) or block2 (ones)- in this expt block2 are LED trials
        b2Ix = zeros(1,length(b2IxRaw));
        for trial = 1:length(b2IxRaw)
            if isnan(b2IxRaw(trial)) == 1;
                b2Ix(trial) = 0;
            elseif b2IxRaw(trial) == 0;
                b2Ix(trial) = 0;
            elseif b2IxRaw(trial) == 1;
                b2Ix (trial) = 1;
            end
        end
        
        
        %only include contrasts present in both blocks
        if ismember(1, b2Ix) == 1; %if block 2 is present, include contrasts from both blocks
            all_cons1 = tGratingContrast(b2Ix==0);
            all_cons2 = tGratingContrast(b2Ix==1);
            cons_both = intersect(all_cons1, all_cons2);
        else %if block 2 isn't present, cons_both is just the block 1 contrasts
            all_cons1 = tGratingContrast(b2Ix==0);
            cons_both = tGratingContrast(b2Ix==0);
        end
        
        % Percent correct vs contrast for Block 1
        %tCons1 = cons_both;
        %cons1 = unique(tCons1); %list of target contrasts
        cons1 = unique(all_cons1);
        ncon1 = length(cons1); %number of target contrasts
        pctCorr1 = zeros(1,ncon1); %empty matrix for percent correct for each target contrast
        for icon = 1:ncon1 %for loop to go through each target
            all_trial_ind = intersect(find(b2Ix==0),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B1 
            nS = length(intersect(find(SIx),all_trial_ind)) %number of correct trials for this target con
            nF = length(intersect(find(MIx),all_trial_ind)) %number of incorrect trials for this target con
            pctCorr1(1,icon) = nS/(nS+nF);
        end
         
        % Percent correct vs contrast for Block 2
        if ismember(1, b2Ix) == 1;
        %tCons2 = cons_both;
        %cons2 = unique(tCons2); %list of target contrasts
        cons2 = unique(all_cons2);
        ncon2 = length(cons2); %number of target contrasts
        pctCorr2 = zeros(1,ncon2); %empty matrix for percent correct for each target contrast
        for icon = 1:ncon2 %for loop to go through each target
            all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons2(icon))); %find all trials for this target con in B1 
            nS = length(intersect(find(SIx),all_trial_ind)); %number of correct trials for this target con
            nF = length(intersect(find(MIx),all_trial_ind)); %number of incorrect trials for this target con
            pctCorr2(1,icon) = nS/(nS+nF);
        end
        end
        
        %{
        mouseText = strjoin(string(miceList)', ', ')
        NSessions = string(length(find(all_data.trialsSinceReset > 0)))
        NTrials = string(length(all_data.trialOutcomeCell))
        %}
        
        % Plot these
        figure(i);
        %sgtitle(strjoin([[mouseText], ", NSessions = ", NSessions, ", NTrials = ", NTrials, ", Threshold Hit Rate = ", thresholdHitRate], ''));
        subplot(3,interval,interval+ifile);
        line1 = plot(cons1(isfinite(pctCorr1)),pctCorr1(isfinite(pctCorr1)),'-ok', 'HandleVisibility','off');
       % set(line1, 'IconDisplayStyle', 'off')
        hold on;
        if ismember(1, b2Ix) == 1;
            line2 = plot(cons2(isfinite(pctCorr2)),pctCorr2(isfinite(pctCorr2)),'-ob', 'HandleVisibility', 'off');
            %set(line2, 'IconDisplayStyle', 'off')
        end
        xlabel('Contrast')
        ylabel('Hit Rate')
        set(gca,'XScale','log')
        ylim([0 1])
        xlim([0 1])
    
        %Error bars
    
        ciPctCorrBounds = zeros(ncon1,2); %another empty matrix to house CI values
        for icon = 1:ncon1;
            all_trial_ind = intersect(find(b2Ix == 0), find(tGratingContrast == cons1(icon)));
            nS = length(intersect(all_trial_ind, find(SIx)));
            nF = length(intersect(all_trial_ind, find(MIx)));
            [~, ciPctCorrBounds(icon,:)] = binofit(nS, (nS+nF));
        end
        ciPctCorr1 = ciPctCorrBounds(:,2) - ciPctCorrBounds(:,1);
        
        %Calculate Block 2 error bars
        if ismember(1, b2Ix) == 1;
            ciPctCorrBounds = zeros(ncon2,2); %another empty matrix to house CI values
            for icon = 1:ncon2;
                all_trial_ind = intersect(find(b2Ix == 1), find(tGratingContrast == cons2(icon)));
                nS = length(intersect(all_trial_ind, find(SIx)));
                nF = length(intersect(all_trial_ind, find(MIx)));
                [~, ciPctCorrBounds(icon,:)] = binofit(nS, (nS+nF));
            end
            ciPctCorr2 = ciPctCorrBounds(:,2) - ciPctCorrBounds(:,1);
        end
        
        % Add error bars to plot
        figure(i);
        subplot(3,interval,interval+ifile);
        hold on;
        err1 = errorbar(cons1, pctCorr1, ciPctCorr1./2, '-ok');
        set(err1, 'DisplayName', 'Iso')
        if ismember(1, b2Ix) == 1;
            err2 = errorbar(cons2, pctCorr2, ciPctCorr2./2, '-ob');
            set(err2, 'DisplayName', 'Cross')
        end
        if ifile == 1;
            Lgnd = legend('show');
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.5;
        end
    
    
    %% Contrast vs reaction time
    
    
        %Get react times
        tRT = celleqel2mat_padded(input.reactTimesMs);
        
        %Empty vectors for block 1 and block 2 average RTs
        AvgRT_1 = zeros(1, ncon1);
        if ismember(1, b2Ix) == 1;
            AvgRT_2 = zeros(1,ncon2);
        end
        
        %Empty vectors for block 1 and block 2 error bars (standard errors)
        SeRT_1 = zeros(1, ncon1);
        if ismember(1, b2Ix) == 1;
            SeRT_2 = zeros(1, ncon2);
        end
        
        for icon = 1:ncon1
            all_trial_ind = intersect(find(b2Ix==0),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B1 
            success_ind = intersect(all_trial_ind, find(SIx));
            RTList = tRT(success_ind);
            AvgRT = mean(RTList);
            stdev = std(RTList);
            se = stdev/sqrt(length(RTList));
            AvgRT_1(:,icon) = AvgRT;
            SeRT_1(:,icon) = se;
        end
        
        
        if ismember(1, b2Ix) == 1;
            for icon = 1:ncon2
                all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B2 
                success_ind = intersect(all_trial_ind, find(SIx));
                RTList = tRT(success_ind);
                AvgRT = mean(RTList);
                stdev = std(RTList);
                se = stdev/sqrt(length(RTList));
                AvgRT_2(:,icon) = AvgRT;
                SeRT_2(:,icon) = se;
            end
        end
        
        
        
        figure(i); 
        subplot(3,interval,2*(interval)+ifile);
        plot(cons1, AvgRT_1, '-ok', 'HandleVisibility','off')
        hold on;
        rt1 = errorbar(cons1, AvgRT_1, SeRT_1, '-ok');
        set(rt1, 'DisplayName', 'Iso')
        if ismember(1, b2Ix) == 1;
            plot(cons2, AvgRT_2, '-ob', 'HandleVisibility','off')
            rt2 = errorbar(cons2, AvgRT_2, SeRT_2, '-ob');
            set(rt2, 'DisplayName', 'Cross')
        end
        xlabel('Contrast')
        ylabel('RT')
        %{
        if i == 1
            Lgnd2 = legend('Show')
            Lgnd.Position(1) = 0.02;
            Lgnd.Position(2) = 0.28;
        end
        %}
    
        
    
    end

%% training mode    
elseif mode == "train"
    dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
    for ifile = 1:length(recentSessions);
       try %To catch errors caused by corrupt files
            load(fullfile(dir,recentSessions(ifile)));
       catch
            disp("One or more files were corrupt and did not load")
            continue
        end 
           
    
        %pull put performance variables
        cons = celleqel2mat_padded(input.tGratingContrast);
        conlist = unique(cons); %list of contrasts presented in the session
        maxCon = max(conlist); %finds the max contrast (easiest contrast to see)
        SIx = strcmp(input.trialOutcomeCell,'success'); %Index of successes (ones)
        MIx = strcmp(input.trialOutcomeCell,'ignore'); %Index of misses (ones)
        FAIx = strcmp(input.trialOutcomeCell,'failure'); %Index of false alarms (ones)
    
    
        %checks hit rate for trials with the highest contrast
        easyTrials = find(cons == maxCon);
        hitOnEasiest = length(intersect(find(SIx),easyTrials));
        missOnEasiest = length(intersect(find(MIx),easyTrials));
        hitRateOnEasiest = hitOnEasiest/(hitOnEasiest + missOnEasiest);
        hitRateSummary(1, ifile) = hitOnEasiest;
        hitRateSummary(2, ifile) = missOnEasiest;
        hitRateSummary(3, ifile) = hitRateOnEasiest;
        
    
        %add the sessions that meet your criteria to a temporary structure
        if requireBlock2 == 1;
            if (block2On == 1) && (hitRateOnEasiest >= thresholdHitRate);
                temp(ifile) = input; 
            else
                continue
            end
        elseif requireBlock2 == 0;
            if hitRateOnEasiest >= thresholdHitRate;
                try %Again, to get around corrupt files
                    temp(ifile) = input; %Structs from imaging sessions are different from rgular structs, and this won't join them
                catch
                    continue
                end 
                    %}
            else
                continue
            end
        end
    
    
    
    %{
    %Concatenate the data in the temporary 1 x n structure 
    % into a single 1 x 1 with all of your selected trials
    all_data = concatenateDataBlocks(temp);
    %}
    
    %% plot online performance
    
    
        % Online performance
        trialNum = linspace(1, length(input.tGratingContrast), length(input.tGratingContrast));
        nTrial = length(trialNum);
    
      
    
        dateText = string(input.saveTime);
        fixedTime = string(input.fixedReqHoldTimeMs);
        randTime = string(input.randReqHoldMaxMs);
        tooFastTime = string(input.tooFastTimeMs);
        parameterText = join(["Fixed =", fixedTime, "Random =", randTime, "TooFast =", tooFastTime]);
    
        figure(i);
        sgtitle(mouse);
        subplot(3,interval,ifile)
        hold on;
        if ifile <= 2
        perf1 = plot(smooth(SIx, ceil(nTrial/10), 'lowess'));
        ylim([0 1]);
       
        perf2 = plot(smooth(SIx, nTrial, 'lowess'));
        set(perf2, 'Color', 'r', ...
                'Linewidth', 3)
        ylim([0 1]);
        
        perf3 = plot(smooth(SIx, 100, 'lowess'));
        set(perf3, 'Color', 'k', ...
            'LineWidth', 2);
        ylim([0 1]);
         
        perf4 = plot(smooth(MIx, 100, 'lowess'));
        set(perf4, 'Color', 'm', ...
                'LineWidth', 2, ...
                'LineStyle', '-.');
        ylim([0 1]);
       
        perf5 = plot(smooth(FAIx, 100, 'lowess'));
        set(perf5, 'Color', [0 0.6 0], ...
                'LineWidth', 1, ...
                'LineStyle', '-.');
        ylim([0 1]);

        xlabel('Trial Number')
        ylabel('Rate')
        fig = gcf;
        fig.Position(3) = fig.Position(3) + 250;
        title(dateText)
        subtitle(parameterText, "FontSize", 8)
        else
        perf1 = plot(smooth(SIx, ceil(nTrial/10), 'lowess'));
        set(perf1, 'DisplayName', 'Hit 1');
        ylim([0 1]);
       
        perf2 = plot(smooth(SIx, nTrial, 'lowess'));
        set(perf2, 'Color', 'r', ...
                'Linewidth', 3, 'DisplayName', 'Hit 2')
        ylim([0 1]);
        
        perf3 = plot(smooth(SIx, 100, 'lowess'));
        set(perf3, 'Color', 'k', ...
            'LineWidth', 2, 'DisplayName', 'Hit 3');
        ylim([0 1]);
         
        perf4 = plot(smooth(MIx, 100, 'lowess'));
        set(perf4, 'Color', 'm', ...
                'LineWidth', 2, ...
                'LineStyle', '-.', 'DisplayName', 'Miss');
        ylim([0 1]);
       
        perf5 = plot(smooth(FAIx, 100, 'lowess'));
        set(perf5, 'Color', [0 0.6 0], ...
                'LineWidth', 1, ...
                'LineStyle', '-.', 'DisplayName', 'FA');
        ylim([0 1]);

        xlabel('Trial Number')
        ylabel('Rate')
        title(dateText)
        subtitle(parameterText, "FontSize", 8)
    
        Lgnd = legend('show');
        Lgnd.Position(1) = 0.02;
        Lgnd.Position(2) = 0.73;
            %legend('Hit', 'Hit', 'Hit(Mid)', 'Miss', 'FA', 'Location','eastoutside')
        end
    
    %%
    %hold time histogram
%{
   holdTimes = celleqel2mat_padded(input.holdTimesMs);
   avgRewardedTime = sum(celleqel2mat_padded(input.tTotalReqHoldTimeMs))/length(input.tTotalReqHoldTimeMs)

   figure(i)
   subplot(3,interval,(interval)+ifile);
   histogram(holdTimes, 20, 'BinLimits',[0,5000])
   hold on
   xlabel('Hold Time')
   ylabel('Frequency')
   xline(avgRewardedTime, '--k', 'LineWidth', 1)
%}

    %% hold time plot
%{
axH = subplot(3,interval,(interval)+ifile);
hold on;

trXLim = [0 nTrial]; %get(gca, 'XLim');
reactV = celleqel2mat_padded(input.reactTimesMs);
holdV = celleqel2mat_padded(input.holdTimesMs);
dnscIx = celleqel2mat_padded(input.tDoNoStimulusChange)==1

    % do smooth on only corrects and earlies, then plot by true trial number
desIx = SIx|FAIx;
xVals = find(desIx);
v1 = smooth(holdV(desIx), 25, 'rloess');
v2 = smooth(holdV(desIx), 250, 'rlowess');
desYIx = SIx;
xYVals = find(desYIx);
vy1 = smooth(reactV(SIx), 25, 'rloess');
vy2 = smooth(reactV(SIx), 250, 'rlowess');
[axH pH1 pH2] = plotyy(1,1,1,1);
set(axH, 'NextPlot', 'add');

% first axes
if ~isempty(v1) && ~isempty(v2)
  hH(1) = plot(axH(1), xVals, v1);
  hH(2) = plot(axH(1), xVals, v2);
  
  c1 = 'k';
  set(hH, 'Color', c1);
  set(hH(2), ...
           'LineWidth', 3);

  set(axH, 'XLim', trXLim, ...
           'YLimMode', 'auto', ...
           'YTickMode', 'auto', ...
           'YTickLabelMode', 'auto', ...
           'YColor', c1);
  yLim = get(axH(1), 'YLim');
  if yLim(2) > 0
    yLim(1) = 0;
    set(axH(1), 'YLim', yLim);
  end
  if ifile == 1
    ylabel('Hold time (ms) - corr+early');
  end
end

% 2nd axes
if ~isempty(vy1) && ~isempty(vy2)
    hyH(1) = plot(axH(2), xYVals, vy1);
    hyH(2) = plot(axH(2), xYVals, vy2);
    if isfield(input, 'tDoNoStimulusChange') && sum(dnscIx)>0,
        dnscCorrIx = dnscIx & SIx;
        dnscCorrIx = dnscCorrIx(dnscIx);
        dnscVals = smooth(reactV(dnscCorrIx), 25, 'rloess');
        plot(axH(2), dnscCorrTrs, dnscVals, 'Color', [0 0.6 0]);
    end
    hyH(2) = plot(axH(2), xYVals, vy2);
    
    c2 = 'b';
    set(hyH, 'Color', c2);
    set(hyH(2), 'LineWidth', 2);
    
    set(axH(2), 'YLim', [0 max(reactV)], ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto', ...
                'YColor', c2)
    if ifile == interval;
        ylabel(axH(2), 'React time (ms) - corr');
    end
end

title('mean react/hold over time');
%}
    %% react time histogram

reactTimes = celleqel2mat_padded(input.reactTimesMs);
tooFastTime = input.tooFastTimeMs
[counts, edges] = histcounts(reactTimes, 60, 'BinLimits', [-2000, 2000], 'Normalization','percentage')
cumulCounts = cumsum(counts)
edgesTrimmed = edges(2:61)

figure(i)
   subplot(3,interval,(interval)+ifile);
   histogram(reactTimes, 60, 'BinLimits',[-2000,2000], 'Normalization', 'percentage')
   hold on
   %cdfplot(reactTimes)
   plot(edgesTrimmed, cumulCounts)
   xlabel('React Time')
   ylabel('Cumulative Percentage')
   ylim([0,105])
   %xline(tooFastTime, '--k', 'LineWidth', 1)
    


%% work over time plot
clear firstRewMs
clear totalRewMs

axH = subplot(3,interval,2*(interval)+ifile);
hold on;
holdStarts = double(cellvect2mat_padded(input.holdStartsMs));
hSDiffsSec = diff(holdStarts)/1000;
smoothType = 'lowess'
trXLim = [0 nTrial];

% make outliers a fixed value
largeIx = hSDiffsSec >= 120;
hSCapped = hSDiffsSec;  
hSCapped(largeIx) = 120;  

xs = 1:length(hSDiffsSec);




if ~isempty(hSDiffsSec) && sum(~isnan(hSDiffsSec)) > 1
  % computations here
  trPerMin = 1./hSDiffsSec.*60;

  nTrs = length(input.juiceTimesMsCell);
  firstSize = repmat(NaN, [1 nTrs]);
  totalSize = repmat(NaN, [1 nTrs]);
  for iT=1:nTrs %CHANGE TO 1
    tJ = input.juiceTimesMsCell{iT};
    if isempty(tJ), tJ = 0; end
    firstRewMs(iT) = tJ(1);
    totalRewMs(iT) = sum(tJ);
  end

  rateNSmooth = 45;
  avgRatePerMin = smooth(trPerMin, rateNSmooth, smoothType);
  avgCorrPerMin = smooth(trPerMin.*SIx(2:end), rateNSmooth*1.5, smoothType);
  avgRewPerMin = smooth(trPerMin.*totalRewMs(2:end), rateNSmooth*1.5, smoothType);
  

[axesH pH1 pH1a]=plotyy(xs,hSCapped, xs, avgRatePerMin);
  set(pH1, 'LineStyle', 'none', ...
           'Marker', 'x');
  set(axesH, 'NextPlot', 'add');


  
  if sum(largeIx) > 0
    pH2 = plot(axesH(1),xs(largeIx),hSCapped(largeIx),'r.');  % outliers
  end
  

  % corrects
  pH1b = plot(axesH(2), xs, avgCorrPerMin, 'k');

  
  % avg rew
  pH1b = plot(axesH(2), xs, avgRewPerMin./10/6, 'k--');
  
  
else
  axesH(1) = gca;
  axesH(2) = NaN;
end

if ifile == 1
    ylabel('trial start time diff (s)');
end
xLim = trXLim;
set(axesH, 'XLim', xLim);
lH = plot(xLim, 20*[1 1], '--k');

set(axesH(1),'YLim', [0 121], ...
             'YTick', [0 40 80 120], ...
             'YTickLabelMode', 'auto', ...
             'YColor', 'k');


if ~isempty(axesH(2))
    if ifile == interval
        ylabel(axesH(2), 'Trials/min; avg rew (ms/sec)');
    end
  yLim2 = get(axesH(2), 'YLim');
  yLim2 = [0 max(yLim2(2), 6)];
  set(axesH(2), 'YLim', yLim2, ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto');

end

  
nDiffs = length(hSDiffsSec);
fN = max(1, nDiffs-5);  % if first 6 trials, start at 1
title(sprintf('Last 6 (sec): %s', mat2str(round(hSDiffsSec(fN:end)))));
end
end

%%

today = datetime('now')
yearval = year(today)
monthval = month(today)
dayval = day(today)

if monthval <=9
    monthstr = join(["0", string(monthval)],'')
else
    monthstr = string(monthval)
end

if dayval <=9
    daystr = join(["0", string(dayval)],'')
else
    daystr = string(dayval)
end

yearchar = char(string(yearval))
yearstr = string(yearchar(3:4))

datestr = join([yearstr, monthstr, daystr], '')
savename = join([datestr, mouse], '_')

fig = figure(i)
fig.WindowState = 'maximized'

if saveImage == 1
    saveas(fig, savename,'jpg' )
end
end

%% optional hold time plot for additional animals

%{

 holdTimes = celleqel2mat_padded(input.holdTimesMs);
 avgRewardedTime = sum(celleqel2mat_padded(input.tTotalReqHoldTimeMs))/length(input.tTotalReqHoldTimeMs)


for i = 1:interval
   figure(5)
   subplot(interval,1,i);
   histogram(holdTimes, 20, 'BinLimits',[0,10000])
   hold on
   xlabel('Hold Time')
   ylabel('Frequency')
   xline(avgRewardedTime, '--k', 'LineWidth', 1)
end

%}


