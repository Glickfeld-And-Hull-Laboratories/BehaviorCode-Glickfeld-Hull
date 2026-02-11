
%%Load Data

% if you try to load sessions with block 1 together with sessions without
% block 1, the whole b2Ix list comes back undefined

clear all
%Import excel spreadsheet with file info
%(note that you may have to update the columns you use as you go using the "Range" parameter)
excel = readtable("home\LindseyW\MATLAB\Behavior_analysis\BehaviorIndices.xlsx", Range="A363:D380", Sheet = 1, ReadVariableNames = true);

selectMice = 1; %Set to zero if you want to load data from all mice, set to 1 if you want to load data from specific mice
miceList = ['i2131']; %If selectMice = 1, this list contains the mice whose data you want to load (separate mice ID's with semicolons)

%Find sessions you want to load
mouse = excel{:,2};
if selectMice == 1;
    goodSessions = find(ismember(mouse, miceList));
    filenames = string(excel{goodSessions,1});
elseif selectMice == 0;
   filenames = excel{:,1};
end

hitRateSummary = zeros(3, length(filenames));

dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
for ifile = 1:length(filenames);
    try %To catch errors caused b corrupt files
        load(fullfile(dir,filenames(ifile)));
    catch
        continue
    end


    %Checks whether block 1 is on
    b2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
    block2On = input.doBlock2;

    cons = celleqel2mat_padded(input.tGratingContrast);
    conlist = unique(cons); %list of contrasts presented in the session
    maxCon = max(conlist); %finds the max contrast (easiest contrast to see)
    SIx = strcmp(input.trialOutcomeCell,'success'); %Index of successes (ones)
    MIx = strcmp(input.trialOutcomeCell,'failure'); %Index of misses (ones)

    %checks hit rate for trials with the highest contrast
    easyTrials = find(cons == maxCon);
    hitOnEasiest = length(intersect(find(SIx),easyTrials));
    missOnEasiest = length(intersect(find(MIx),easyTrials));
    hitRateOnEasiest = hitOnEasiest/(hitOnEasiest + missOnEasiest);
    hitRateSummary(1, ifile) = hitOnEasiest;
    hitRateSummary(2, ifile) = missOnEasiest;
    hitRateSummary(3, ifile) = hitRateOnEasiest;
    

    
    thresholdHitRate = 0; %set threshold for hit rate on easiest trials, sessions with hit rate below threshold won't be included
    requireBlock2 = 0; %set to 1 if you only want to inlude sessions with block 2 enabled, otherwise set to 0 to include all sessions

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
                temp(ifile) = input;
            catch
                continue
            end
        else
            continue
        end
    end
end

%Concatenate the data in the temporary 1 x n structure 
% into a single 1 x 1 with all of your selected trials
all_data = concatenateDataBlocks(temp);

%%  Contrast vs.Hit Rate

%clearvars -except all_data
close all

% Extract trial outcomes
tGratingContrast = celleqel2mat_padded(all_data.tGratingContrast); %contrast for each trial
SIx = strcmp(all_data.trialOutcomeCell,'success'); %Index of successes (ones)
MIx = strcmp(all_data.trialOutcomeCell,'ignore'); %Index of misses (ones) (wouldn't misses be coded as 'ignore'?)
b2IxRaw = celleqel2mat_padded(all_data.tBlock2TrialNumber); %IDs trial as block1 (zeros) or block2 (ones)- in this expt block2 are LED trials
for trial = 1:length(b2IxRaw)
    if isnan(b2IxRaw) == 1;
        b2Ix(trial) = 0;
    elseif b2IxRaw == 0;
        b2Ix (trial) = 0;
    elseif b2IxRaw == 1;
        b2Ix (trial) = 1;
    end
end


%only include contrasts present in both blocks
if ismember(1, b2Ix) == 1; %if block 2 is present, include contrasts from both blocks
    all_cons1 = tGratingContrast(b2Ix==0);
    all_cons2 = tGratingContrast(b2Ix==1);
    cons_both = intersect(all_cons1, all_cons2);
else %if block 2 isn't present, cons_both is just the block 1 contrasts
    cons_both = tGratingContrast(b2Ix==0);
end

% Percent correct vs contrast for Block 1
tCons1 = cons_both;
cons1 = unique(tCons1); %list of target contrasts
ncon1 = length(cons1); %number of target contrasts
pctCorr1 = zeros(1,ncon1); %empty matrix for percent correct for each target contrast
for icon = 1:ncon1 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==0),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B1 
    nS = length(intersect(find(SIx),all_trial_ind)); %number of correct trials for this target con
    nF = length(intersect(find(MIx),all_trial_ind)); %number of incorrect trials for this target con
    pctCorr1(1,icon) = nS/(nS+nF);
end
 
% Percent correct vs contrast for Block 2
if ismember(1, b2Ix) == 1;
tCons2 = cons_both;
cons2 = unique(tCons2); %list of target contrasts
ncon2 = length(cons2); %number of target contrasts
pctCorr2 = zeros(1,ncon2); %empty matrix for percent correct for each target contrast
for icon = 1:ncon2 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons2(icon))); %find all trials for this target con in B1 
    nS = length(intersect(find(SIx),all_trial_ind)); %number of correct trials for this target con
    nF = length(intersect(find(MIx),all_trial_ind)); %number of incorrect trials for this target con
    pctCorr2(1,icon) = nS/(nS+nF);
end
end

mouseText = strjoin(string(miceList)', ', ')
NSessions = string(length(find(all_data.trialsSinceReset > 0)))
NTrials = string(length(all_data.trialOutcomeCell))

% Plot these
figure(1);
sgtitle(strjoin([[mouseText], ", NSessions = ", NSessions, ", NTrials = ", NTrials, ", Threshold Hit Rate = ", thresholdHitRate], ''));
subplot(2,2,1);
plot(cons1(isfinite(pctCorr1)),pctCorr1(isfinite(pctCorr1)),'ok');
hold on;
if ismember(1, b2Ix) == 1;
    plot(cons2(isfinite(pctCorr2)),pctCorr2(isfinite(pctCorr2)),'ob');
end
xlabel('Contrast')
ylabel('Hit Rate')
set(gca,'XScale','log')
ylim([0 1])
xlim([0 1])




% Fit cumulative Weibull to curve
%wbl_1 = fitdist(pctCorr1', 'wbl', 'by', cons1');



%Calculate Block 1 error bars
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
figure(1);
subplot(2,2,1);
hold on;
errorbar(cons1, pctCorr1, ciPctCorr1./2, '-ok');
if ismember(1, b2Ix) == 1;
errorbar(cons2, pctCorr2, ciPctCorr2./2, '-ob');
end
legend('Iso', 'Cross', 'Location', 'southeast')






%% Cumulative hits, misses, and false alarms over trial number
%{
%Trial number vector
trialNum = linspace(1, length(tGratingContrast), length(tGratingContrast));

%Trial outcome vectors
tHit = SIx;
tMiss = MIx;
tFA = strcmp(all_data.trialOutcomeCell,'failure');

totalHit = zeros(1, length(tGratingContrast));
totalMiss = zeros(1, length(tGratingContrast));
totalFA = zeros(1, length(tGratingContrast));

holdTimes = all_data.holdTimesMs;
reqTimes = all_data.tTotalReqHoldTimeMs;

for itrial = 1:length(tGratingContrast)
    totalHit = cumsum(tHit);
    totalMiss = cumsum(tMiss);
    totalFA = cumsum(tFA);
end

%Plot cumulative vectors over trial number
figure(2);
hold on;
xlabel('Trial number');
ylabel('Cumulative count');
plot(trialNum, totalHit, '-k');
plot(trialNum, totalMiss, '-b');
plot(trialNum, totalFA, '-g');
legend('Hits', 'Misses', 'False Alarms', 'Location', 'northeastoutside');
%}

%% Pct Correct vs. trial number
%{
trialNum = linspace(1, length(input.tGratingContrast), length(input.tGratingContrast));
nTrial = length(trialNum)

inputHit = strcmp(input.trialOutcomeCell,'success');
inputMiss = strcmp(input.trialOutcomeCell,'ignore');
inputFA = strcmp(input.trialOutcomeCell,'failure');

%if input.doBlock2;
   % input




figure(3);
hold on;
perf1 = plot(smooth(inputHit, ceil(nTrial/10), 'lowess'));
perf2 = plot(smooth(inputHit, nTrial, 'lowess'));
set(perf2, 'Color', 'r', ...
            'Linewidth', 3)
perf3 = plot(smooth(inputHit, 100, 'lowess'));
set(perf3, 'Color', 'k', ...
        'LineWidth', 2);
perf4 = plot(smooth(inputMiss, 100, 'lowess'));
set(perf4, 'Color', 'm', ...
            'LineWidth', 2, ...
            'LineStyle', '-.');



%for itrial = 1:length(trialNum)
%}

%% Contrast vs. FA rate

% Extract trial outcomes
FAIx = strcmp(all_data.trialOutcomeCell,'failure');


% Percent correct vs contrast for Block 1
tCons1 = cons_both;
cons1 = unique(tCons1); %list of target contrasts
ncon1 = length(cons1); %number of target contrasts
pctFA_1 = zeros(1,ncon1); %empty matrix for percent correct for each target contrast
for icon = 1:ncon1 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==0),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B1 
    nFA = length(intersect(find(FAIx),all_trial_ind)); %number of correct trials for this target con
    n_notFA = length(intersect(find(MIx),all_trial_ind)) + length(intersect(find(SIx),all_trial_ind)); %number of incorrect trials for this target con
    pctFA_1(1,icon) = nFA/(nFA+n_notFA);
end
 
% Percent correct vs contrast for Block 2
if ismember(1, b2Ix) == 1;
tCons2 = cons_both;
cons2 = unique(tCons2); %list of target contrasts
ncon2 = length(cons2); %number of target contrasts
pctFA_2 = zeros(1,ncon2); %empty matrix for percent correct for each target contrast
for icon = 1:ncon2 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons2(icon))); %find all trials for this target con in B1 
    nFA = length(intersect(find(FAIx),all_trial_ind)); %number of correct trials for this target con
    n_notFA = length(intersect(find(MIx),all_trial_ind)) + length(intersect(find(SIx),all_trial_ind)); %number of incorrect trials for this target con
    pctFA_2(1,icon) = nFA/(nFA+n_notFA);
end
end

% Plot these
figure(1);
subplot(2,2,2)
plot(cons1(isfinite(pctFA_1)),pctFA_1(isfinite(pctFA_1)),'-ok');
hold on
if ismember(1, b2Ix) == 1;
plot(cons2(isfinite(pctFA_2)),pctFA_2(isfinite(pctFA_2)),'-ob');
end
xlabel('Contrast')
ylabel('Percent False Alarm')
set(gca,'XScale','log')
ylim([0 1])
xlim([0 1])
legend('Iso', 'Cross')

%% Contrast vs. Miss rate
%{
% Percent correct vs contrast for Block 1
tCons1 = cons_both;
cons1 = unique(tCons1); %list of target contrasts
ncon1 = length(cons1); %number of target contrasts
pctMiss_1 = zeros(1,ncon1); %empty matrix for percent correct for each target contrast
for icon = 1:ncon1 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==0),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B1 
    nMiss = length(intersect(find(MIx),all_trial_ind)); %number of correct trials for this target con
    n_notMiss = length(intersect(find(FAIx),all_trial_ind)) + length(intersect(find(SIx),all_trial_ind)); %number of incorrect trials for this target con
    pctMiss_1(1,icon) = nMiss/(nMiss+n_notMiss);
end
 
% Percent correct vs contrast for Block 2
tCons2 = cons_both;
cons2 = unique(tCons2); %list of target contrasts
ncon2 = length(cons2); %number of target contrasts
pctMiss_2 = zeros(1,ncon2); %empty matrix for percent correct for each target contrast
for icon = 1:ncon2 %for loop to go through each target
    all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons2(icon))); %find all trials for this target con in B1 
    nMiss = length(intersect(find(MIx),all_trial_ind)); %number of correct trials for this target con
    n_notMiss = length(intersect(find(FAIx),all_trial_ind)) + length(intersect(find(SIx),all_trial_ind)); %number of incorrect trials for this target con
    pctMiss_2(1,icon) = nMiss/(nMiss+n_notMiss);
end

% Plot these
figure(1); 
subplot(2,3,3)
plot(cons1(isfinite(pctMiss_1)),pctMiss_1(isfinite(pctMiss_1)),'-ok');
hold on;
plot(cons2(isfinite(pctMiss_2)),pctMiss_2(isfinite(pctMiss_2)),'-ob');
xlabel('Contrast')
ylabel('Percent Miss')
set(gca,'XScale','log')
ylim([0 1])
xlim([0 1])
legend('Iso', 'Cross')
%}

%% Contrast vs. Reaction time

%Get react times
tRT = celleqel2mat_padded(all_data.reactTimesMs);

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
    success_ind = intersect(all_trial_ind, find(SIx))
    RTList = tRT(success_ind);
    AvgRT = mean(RTList);
    stdev = std(RTList);
    se = stdev/sqrt(length(RTList))
    AvgRT_1(:,icon) = AvgRT;
    SeRT_1(:,icon) = se
end

if ismember(1, b2Ix) == 1;
for icon = 1:ncon2
    all_trial_ind = intersect(find(b2Ix==1),find(tGratingContrast==cons1(icon))); %find all trials for this target con in B2 
    success_ind = intersect(all_trial_ind, find(SIx))
    RTList = tRT(success_ind);
    AvgRT = mean(RTList);
    stdev = std(RTList);
    se = stdev/sqrt(length(RTList))
    AvgRT_2(:,icon) = AvgRT;
    SeRT_2(:,icon) = se
end
end



figure(1); 
subplot(2,2,3);
plot(cons1, AvgRT_1, 'ok')
hold on;
errorbar(cons1, AvgRT_1, SeRT_1, '-ok');
if ismember(1, b2Ix) == 1;
    plot(cons2, AvgRT_2, 'ob')
    errorbar(cons2, AvgRT_2, SeRT_2, '-ob');
end
xlabel('Contrast')
ylabel('RT')
legend('Iso', 'Cross')




%% Reaction time t-test

%Calculate overall averages
if ismember(1, b2Ix) == 1;
    OverallAvgRT_1 = mean(AvgRT_1)
    OverallAvgRT_2 = mean(AvgRT_2)
    OverallAvgRT_both = [OverallAvgRT_1,OverallAvgRT_2]
else
     OverallAvgRT_both = mean(AvgRT_1)
end
%Set x -axis for bar graph
x = categorical({'Iso','Cross'})

%Calculate overall error bars (standard errors)
success_ind_1 = intersect(find(b2Ix==0),find(SIx))
RTs_1 = tRT(success_ind_1)
OverallSeRT_1 = std(RTs_1)/sqrt(length(RTs_1))
OverallSeRT_both = OverallSeRT_1

if ismember(1, b2Ix) == 1;
    success_ind_2 = intersect(find(b2Ix==1),find(SIx))
    RTs_2 = tRT(success_ind_2)
    OverallSeRT_2 = std(RTs_2)/sqrt(length(RTs_2))
    OverallSeRT_both = [OverallSeRT_1, OverallSeRT_2]
%T-test to see if overall error bars are different
    [h,p,ci,stats] = ttest2(RTs_1, RTs_2)
    pvalText = strjoin(["p =", string(round(p,3))])
end

figure(1);
subplot(2,2,4);
b = bar(x,OverallAvgRT_both)
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 0];
b.CData(2,:) = [0 0 1];
hold on
errorbar(x ,OverallAvgRT_both, OverallSeRT_both, 'or')
ylim([300 400])
xlabel('Stimulus Type')
ylabel('RT')
title(pvalText)

%% 







