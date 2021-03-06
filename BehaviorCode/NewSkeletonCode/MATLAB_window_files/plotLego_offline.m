function input = plotLego(input)

% essential consts
name = 'Lego';
cs = lego_constants;
spSz = {4,3};

smoothType = 'lowess';
figNum = 1;

%% other consts
yColor = [1 1 0]*0.5;

%% draw figure
figH = figure(figNum);
set(figH, 'ToolBar', 'none');
clf;
%set(figH, 'Visible', 'off');

switch hostname
 case 'test-rig'
  figPos = [1111 208 806 923];
 otherwise
  figPos = [700 140 780 770];
end
set(figH, 'Position', figPos);

%% set up arrays
if ~isfield(input, 'changedStr')
  input.changedStr = {};
end

input.savedDataName = sprintf('~/Documents/MWorks/Data/data-i%03d-%s.mat', ...
                              input.subjectNum, input.saveTime);

% some data processing
nTrial = length(input.trialOutcomeCell);
block2Ix = cell2mat_padded(input.tBlock2TrialNumber);

correctIx = strcmp(input.trialOutcomeCell, 'success');
incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
missIx = ignoreIx;
nCorr = sum(correctIx);
nInc = sum(incorrectIx);
nIg = sum(ignoreIx);
nTrials = length(input.trialOutcomeCell);
decisionV = celleqel2mat_padded(input.tDecisionTimeMs);
trPer80Str = regexprep(['[' num2str(input.trPer80V, '%3d') ']'], '[ *', '[');
leftTrPer80Str = ['[' num2str([input.leftTrPer80Level1 input.leftTrPer80Level2 input.leftTrPer80Level3 input.leftTrPer80Level4 input.leftTrPer80Level5 input.leftTrPer80Level6 input.leftTrPer80Level7 input.leftTrPer80Level8]) ']' ];
leftTrialIx = cell2mat_padded(input.tLeftTrial);
leftTrNs = find(leftTrialIx);
rightTrialIx = ~leftTrialIx;
rightTrNs = find(rightTrialIx);
leftTrialIx = logical(leftTrialIx);
nLeft = sum(leftTrialIx);
nRight = sum(rightTrialIx);
block2Ix = cell2mat_padded(input.tBlock2TrialNumber);
noGoIx = celleqel2mat_padded(input.isNoGo);

leftOutcomes = input.trialOutcomeCell(leftTrialIx);
leftCorr = strcmp(leftOutcomes, 'success');
nLeftCorr = sum(leftCorr);
leftIgn = strcmp(leftOutcomes, 'ignore');
nLeftIgn = sum(leftIgn);
leftInc = strcmp(leftOutcomes, 'incorrect');
nLeftInc = sum(leftInc);

rightOutcomes = input.trialOutcomeCell(~leftTrialIx);
rightCorr = strcmp(rightOutcomes, 'success');
nRightCorr = sum(rightCorr);
rightIgn = strcmp(rightOutcomes, 'ignore');
nRightIgn = sum(rightIgn);
rightInc = strcmp(rightOutcomes, 'incorrect');
nRightInc = sum(rightInc);

% Left bias indexing
left_correct_ind = find((double(leftTrialIx)+double(correctIx))==2);
left_incorrect_ind = intersect(find(leftTrialIx==0), find(incorrectIx==1));
tLeftResponse = zeros(size(correctIx));
tLeftResponse([left_correct_ind left_incorrect_ind]) = 1;
tLeftResponse(find(ignoreIx)) = NaN;
input.tLeftResponse = tLeftResponse;

%Right decision time indexing
right_correct_ind = find((double(rightTrialIx)+double(correctIx))==2);


%  TODO:  this should be automated in a loop, or better should not have to convert from cell at all
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;

% stimulus str            
stimStr = strcat(mat2str(input.gratingMaxDiameterDeg), ' deg, ', ...
    num2str(input.gratingSpatialFreqCPD), ' cpd, ', mat2str(input.gratingEccentricityDeg), ' deg');

if input.gratingSpeedDPS > 0,
    stimStr = strcat(stimStr, ', ', mat2str(input.gratingSpeedDPS), ' dps, ', ...
        mat2str(input.gratingStartingPhaseDeg), ' deg');
end

% leftStr formatting
if input.doLeftSeparateOdds == 0
  leftStr = sprintf('leftTrPer80: Matched');
else
  leftStr= sprintf('leftTrPer80: %s \n', leftTrPer80Str);
end
% blockStr formatting
if input.doBlocks==0,
  blockStr = sprintf('Blocks: off, probLeft: %1.2f, probSwitch: %1.2f \n', double(input.tStimProbAvgLeft{nTrials}), double(input.stimProbAvgSwitch));

elseif input.doBlocks==1,
  blockStr = sprintf('Blocks: on, L: %2.0f, R:%2.0f \n', double(input.blockLeftTrs), double(input.blockRightTrs));
end

if input.doBlock2==1,
  block2Str = 'Block2: on \n';
  if input.block2DoTrialLaser==1,
    block2Str = sprintf('Block2: on, Trial Laser Power:%5.2fmW \n', double(input.block2TrialLaserPowerMw));
  end
elseif input.doBlock2==0,
  block2Str = 'Block2: off \n';
end

if isfield(input, 'doMatchToTarget')
  if input.doMatchToTarget==1,
    MTTstr = 'On';
  else
    MTTstr = 'Off';
  end
else
    MTTstr = 'Off';
end

if input.doConsecCorrectReward==1,
  rewStr = strcat(mat2str(input.consecCorrRewardInterval/1000), 'x nCorr');
else
  rewStr = mat2str(input.rewardTimeUs/1000);
end
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf(1,'%d %d\n', sum(block2Tr1Ix&~earlyIx), sum(block2Tr2Ix&~earlyIx));  %% debug # of block2 trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Values

        axH = subplot(spSz{:}, 1);						% default axes are 0 to 1
        
	set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

        numTrials = nTrial;
        
        set(gcf, 'Visible', 'off'); % hide figure during text
                                    % drawing - kludge

        text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 16);

        elMin = round((now - datenum(input.startDateVec)) * 24*60);
        startStr = datestr(input.startDateVec, 'HH:MM');
        text(0.00, 1.05, {'Subject:', 'Start time + elapsed:', 'Reward vol:'});
	text(0.60, 1.05, ...
             { sprintf('%2d', input.subjectNum), ...
               sprintf('%s + %2dm', ...
                       startStr, elMin), ...
               sprintf('%.1f s     \t(%g ms)', ...
                       nansum(juiceTimesMsV./1000), ...
                       chop(nanmean(juiceTimesMsV),2)), ...
             });
         
	t2H(1) = text(0.00, 0.8, {'    ', 'Trials:', 'Correct:', 'Incorrect:', 'Missed:', 'Decision (ms):'});
	t2H(2) = text(0.4, 0.8, {'Left', sprintf('%d', nLeft), sprintf('%d', nLeftCorr), ...
				sprintf('%d', nLeftInc), sprintf('%d', nLeftIgn), ...
                sprintf('%5.0f', median(decisionV(left_correct_ind)))});
	t2H(3) = text(0.65, 0.8, {'Right', sprintf('%d', nRight), sprintf('%d', nRightCorr), ...
				sprintf('%d', nRightInc), sprintf('%d', nRightIgn), ...
                sprintf('%5.0f', median(decisionV(right_correct_ind)))});
    t2H(4) = text(0.9, 0.8, {'Total', sprintf('%d', nTrial), sprintf('%d', nCorr), ...
                sprintf('%d', nInc), sprintf('%d', nIg), ...
                sprintf('%5.0f', median(decisionV(correctIx)))});
    t2H(5) = text(1.15, 0.8, {'  ', '  ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
				sprintf('%.0f%%', nInc / numTrials * 100.0), ...
				sprintf('%.0f%%', nIg / numTrials * 100.0)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');


%         tStr = sprintf(['Decision median:\n', ...
%                         '   %5.1f ms\n'], ...
%                        median(decisionV(correctIx)));
%                
%         text(0.8, 0.8, tStr, ...
%              'VerticalAlignment', 'top', ...
%              'HorizontalAlignment', 'left');
         
if isfield(input, 'gratingTargetEccentricity')
  input.leftDecisionThreshold = (input.gratingEccentricityDeg-input.gratingTargetEccentricity)./ input.feedbackMotionSensitivity;
  input.rightDecisionThreshold = (input.gratingEccentricityDeg-input.gratingTargetEccentricity)./ input.feedbackMotionSensitivity;
end


        tStr = sprintf( ['Decision Time: \t%5.2f s;   ITI %d ms \n', ...
                         'Stim On Time: \t%5d ms \n', ...
                         'Thresholds (L,R):\t%4.0f, %4.0f pulses\n', ...
                         'Too-Fast Time: %d ms;  Reward: %s ms \n', ...
                         'Timeouts (ign,inc):\t%4.1f, %4.1f s\n', ...
                         'trPer80: %s\n', ...
                         leftStr, ...
                         'Stim: %s \n', ...
                         blockStr, ...
                         block2Str, ...
                         'Motion Sensitivity: \t%2.3f \n'], ...
                        input.reactionTimeMs/1000, ...
                        input.itiTimeMs, ...
                        input.stimOnTimeMs, ...
                        input.leftDecisionThreshold, ...
                        input.rightDecisionThreshold, ...
                        input.tooFastTimeMs, ...
                        rewStr, ...
                        input.ignoreTimeoutMs/1000, ...
                        input.incorrectTimeoutMs/1000, ...
                        trPer80Str, ...
                        stimStr, ...
                        input.feedbackMotionSensitivity);

        text(0.0, 0.3, tStr, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

        set(gcf, 'Visible', 'on');
        
%%%%%%%%%%%%%%%

%% 2 - smoothed perf curve
axH = subplot(spSz{:},2);
hold on;
lH2 = plot(smooth(double(correctIx), nTrials/5, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(smooth(double(ignoreIx), nTrials/5, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2);

lH4 = plot(smooth(double(incorrectIx), nTrials/5, smoothType));
  set(lH4, 'Color', 'g', ...
           'LineWidth', 2);
  anystack(lH4, 'bottom');

  anystack(lH3, 'bottom');


title('Outcome Plot');
ylabel('Occurrence (%)');
set(gca, 'YLim', [0 1]);

xlabel('Trials');
trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);

%%%%%%%%%%%%%%%%
  
%% 4 - List changed text


axH = subplot(spSz{:}, 4);
hold on;
set(axH, 'Visible', 'off');
    set(axH, 'OuterPosition', [0.02 0.47, 0.25, 0.2])
if nTrial > 1 
  % recompute text list from all changes
  desIx = ~cellfun(@isempty, input.constChangedOnTrial);
  desIx(1) = false; % skip the first trial
  desTrNs = find(desIx);
  
  changedStrOut = {};
  for iT = 1:length(desTrNs)
    tDesN = desTrNs(iT);
    tChangedList = input.constChangedOnTrial{tDesN};
    [nRows nCols] = size(tChangedList);
    
    trPer80Changed = false;
    for iR = 1:nRows
      tVarName = tChangedList{iR,1};
      tChangedTo = tChangedList{iR,2};
      % special case odds changes
      if strfind(tVarName, 'trPer80Level')
        trPer80Changed = true;  % and iterate; summarize below
      else
          if ischar(tChangedTo)
              tChangedToStr = sprintf('''%s''', tChangedTo);
          else
              tChangedToStr = sprintf('%g', tChangedTo);
          end
          changedStrOut{end+1} = sprintf('Tr %3d - %s: -> %s', ...
                                         tDesN, tVarName, ...
                                         tChangedToStr);
      end
    end
    if trPer80Changed && tDesN > 1
      tStr = mat2str(input.trPer80History{tDesN});
      changedStrOut{end+1} = sprintf('Tr %3d - trPer80LevelN: -> %s', ...
                                     tDesN, tStr);
    end   
  end
  input.changedStr = changedStrOut;

  text(0,1,input.changedStr, ...
       'HorizontalAlignment', 'left', ...
       'VerticalAlignment', 'top');

%%%%%%%%%%

%% 5 - bias plot
axH = subplot(spSz{:},5);
hold on;
if nTrial > 100,
    amtSmooth = round(nTrial*0.10);
else
    amtSmooth = 10;
end

lh1 = plot(smooth(double(leftTrialIx), amtSmooth, smoothType));
set(lh1, 'Color', 'k', 'LineWidth', 2);

lh2 = plot(smooth(double(tLeftResponse), amtSmooth, smoothType));
set(lh2,'Color', 'r', 'LineWidth', 2);

lh3 = refline(0,0.5);
set(lh3, 'LineStyle', '--');
    
title('Bias Plot: Red = Responses, Black = Presentations');
ylabel('Left (%)');
set(gca, 'YLim', [0 1], ...
         'YTick', [0,0.25,0.5,0.75,1]);

xlabel('Trials');
trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%% 3 - Left and Right Outcomes

axH = subplot(spSz{:},3);
hold on;
if sum(leftTrialIx)>0,
  lH2 = plot(leftTrNs, smooth(double(leftCorr), sum(leftTrNs)/5, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(leftTrNs, smooth(double(leftIgn), sum(leftTrNs)/5, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2);

lH4 = plot(leftTrNs, smooth(double(leftInc), sum(leftTrNs)/5, smoothType));
  set(lH4, 'Color', 'g', ...
           'LineWidth', 2);
  anystack(lH4, 'bottom');

  anystack(lH3, 'bottom');
end
if sum(~leftTrialIx)>0,
axH = subplot(spSz{:},3);
hold on;
lH2 = plot(rightTrNs, smooth(double(rightCorr), sum(rightTrNs)/5, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2, 'LineStyle', '--');

lH3 = plot(rightTrNs, smooth(double(rightIgn), sum(rightTrNs)/5, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2, 'LineStyle', '--');

lH4 = plot(rightTrNs, smooth(double(rightInc), sum(rightTrNs)/5, smoothType));
  set(lH4, 'Color', 'g', ...
           'LineWidth', 2, 'LineStyle', '--');
  anystack(lH4, 'bottom');

  anystack(lH3, 'bottom');
  end

title('Outcome Plot ( Left=|, Right=:)');
ylabel('Occurrence (%)');
set(gca, 'YLim', [0 1]);

xlabel('Trials');
trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);

%%%%%%%%%%%%%%%

% 6 - decision time over time
axH = subplot(spSz{:}, 6);
hold on;

correctIx = correctIx & ~noGoIx;
incorrectIx = incorrectIx & ~noGoIx;

noMissIx = correctIx|incorrectIx;
xYVals = find(noMissIx);
decTimes = cell2mat_padded(input.tDecisionTimeMs);
noMissDecTimesMs = double(decTimes(noMissIx));
vy = smooth(noMissDecTimesMs, 25, 'rloess');
decMax = input.reactionTimeMs/1000;

if ~isempty(vy)
    if input.stationaryPeriodMs>0
        startTimes = double(cell2mat(input.tThisTrialStartTimeMs))./1000;
        stimOnMs = double(cell2mat(input.stimTimestampMs))./1000;
        itiTime = (input.itiTimeMs +input.stationaryPeriodMs)./1000;
        delayTime = stimOnMs-startTimes-double(itiTime);
        [axH h1 h2] = plotyy(xYVals, vy,1:nTrial, delayTime);
        set(h1, 'Color', 'b');
        set(h1, 'LineWidth', 2);
        set(axH(1), 'YLim', [0 decMax]);
        ylabel(axH(1), 'Decision time');
        
        set(h2, 'LineStyle', 'none', ...
           'Marker', 'x', 'color', 'k');
        ylabel(axH(2), 'Stationary period (s)');
    else
        hyH = plot(axH, xYVals, vy);
        c2 = 'b';
        set(hyH, 'Color', c2);
        set(hyH, 'LineWidth', 2);

        set(axH, 'YLim', [0 decMax]);
        ylabel(axH, 'Decision time');
    end
    axis tight;
    xlim([1 nTrial]);
end

title('Decision Median over Time');
xlabel('Trials');


%% 7 - Constrast Difference x Correct reaction times plot
axH = subplot(spSz{:},7);
hold on;

if nCorr>0  && input.doTestRobot==0,
  if input.gratingContrastDiffSPO<10
    contDiffV = chop(cell2mat_padded(input.tGratingContrast) ./ cell2mat_padded(input.dGratingContrast),2);
  else
    contDiffV = chop(cell2mat_padded(input.tGratingContrast) - cell2mat_padded(input.dGratingContrast),2);
  end
    corrDiffV = contDiffV(correctIx);
    uqDiff = unique(corrDiffV);
    nLevels = length(uqDiff);
    corrDiffCell = cell(1,nLevels);
    decTimes = cell2mat_padded(input.tDecisionTimeMs);
    corrDecTimes = decTimes(correctIx);

    for ii = 1:length(uqDiff),
      val = uqDiff(ii);
      tix = corrDiffV==val;
      corrDiffCell{ii} = nanmean(corrDecTimes(tix));
    end
    
    if nLeft>0,
      corrDiffVL = contDiffV(correctIx & leftTrialIx);
      uqDiffL = unique(corrDiffVL);
      nLevelsL = length(uqDiffL);
      corrDiffCellL = cell(1,nLevelsL);
      corrDecTimesL = decTimes(correctIx & leftTrialIx);
      for ii = 1:length(uqDiffL),
        val = uqDiffL(ii);
        tix = corrDiffVL==val;
        corrDiffCellL{ii} = nanmean(corrDecTimesL(tix));
      end
    end
    if nRight>0,
      corrDiffVR = contDiffV(correctIx & ~leftTrialIx);
      uqDiffR = unique(corrDiffVR);
      nLevelsR = length(uqDiffR);
      corrDiffCellR = cell(1,nLevelsR);
      corrDecTimesR = decTimes(correctIx & ~leftTrialIx);
      for ii = 1:length(uqDiffR),
        val = uqDiffR(ii);
        tix = corrDiffVR==val;
        corrDiffCellR{ii} = nanmean(corrDecTimesR(tix));
      end
    end


    lev = double(find(input.trPer80V>0)-1);
    if isfield(input, 'dGratingContrastDiff')
      possDiffV = double(input.gratingMaxDiff) ./ (2 .^ (lev./double(input.gratingContrastDiffSPO)))+1;

      minX = min(possDiffV,[],2);
      maxX = max(possDiffV,[],2);
    else
      possDiffV = input.gratingMaxDiff ./ (2 .^ (lev./input.gratingContrastSPO));
      minX = min(possDiffV,[],2);
      maxX = max(possDiffV,[],2);
    end

    minD = min(cell2mat(corrDiffCell));
    maxD = max(cell2mat(corrDiffCell));
    
    hold on
    xLimm = [minX maxX];
    uqDiff = uqDiff;

    if input.gratingContrastDiffSPO<10
      sc_num = 1;
    else
      sc_num = 100;
    end 
    
    pH = plot((sc_num * uqDiff), cell2mat(corrDiffCell));
    if nLeft>0,
      pH1 = plot((sc_num * uqDiffL), cell2mat(corrDiffCellL));
    end
    if nRight>0,
      pH2 = plot((sc_num * uqDiffR), cell2mat(corrDiffCellR));
    end
    if ~(xLimm(1)==0),
      xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
    else
      xL1 = [0 ceil(log10(xLimm(2)))];
    end
    xTickL = 10.^(xL1(1):1:xL1(2));
    xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
    xTLabelL = cellstr(num2str(xTickL(:)));
    

    minX = 0;
    maxX = 100;
    xLimm = [minX maxX];
    if ~(xLimm(1)==0),
      xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
    else
      xL1 = [0 ceil(log10(xLimm(2)))];
    end
    set(gca, 'XLim', [minX maxX], ...
           'XScale', 'log', ...
           'XGrid', 'on');
       set(gca, 'XTick', xTickL);
       set(gca, 'XTickLabel', xTLabelL);


    set(pH, ...
       'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9, ...
        'Color', 'g');
    if nLeft>0,
      set(pH1, ...
          'LineWidth', 1.5, ...
          'Marker', '.', ...
          'MarkerSize', 9, ...
          'Color', 'b')
    end
    if nRight>0,
      set(pH2, ...
          'LineWidth', 1.5, ...
          'Marker', '.', ...
          'MarkerSize', 9, ...
          'Color', [0.8 0.8 0]);
    end

     axis tight  
       ylabel('Decision Time (ms)')
       if input.gratingContrastDiffSPO<10
         xlabel('Contrast Difference (R/L)')
       else
         xlabel('Contrast Difference (R-L)')
      end
       title('Decision Time by Contrast Difference')
end
%%%%%%%%%%%%%%%%%

%% 8 - Contrast Difference x %Correct

axH = subplot(spSz{:},8);
hold on;

if nCorr>0 && input.doTestRobot==0,
    plotTrs = contDiffV(correctIx|incorrectIx);
    uqDiff = unique(plotTrs);
    nLevels = length(uqDiff);
    percentCell = cell(1,nLevels);
    for jj=1:nLevels
        val = uqDiff(jj);
        valIx = contDiffV==val;
        totalNTrialsValB1 = sum(contDiffV(valIx&(correctIx|incorrectIx)));
        corrNTrialsVal = sum(contDiffV(valIx&correctIx));
        percentCell{jj} = corrNTrialsVal/totalNTrialsValB1;
    end
    pH = plot((sc_num * uqDiff), cell2mat(percentCell));
    set(pH, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    set(gca, 'YLim', [0 1], ...
           'XScale', 'log', ...
           'XGrid', 'on');
       set(gca, 'XTick', xTickL);
       set(gca, 'XTickLabel', xTLabelL);
       ylabel('Correct (%)')
       if minX==maxX
         xlim([minX-100 maxX+100])
       elseif minX<maxX
         xlim([minX maxX])
       end
       if input.gratingContrastDiffSPO<10
         xlabel('Contrast Difference (R/L)')
       else
         xlabel('Contrast Difference (R-L)')
      end
       title('Percent Correct by Contrast Difference')
end

%%%%%%%%%%%%%%%%
%% 9 - cdf of decision times

axH = subplot(spSz{:},9);
hold on;

decisionMax = ceil(input.reactionTimeMs/1000)*1000;
if decisionMax <= 6000,
    decisionInterval = 1000;
else
    decisionInterval = 2000;
end
if decisionMax <= 10000,
    decMax = decisionMax;
else
    decMax = 10000;
end

if sum(block2Ix)== 0
  pH = cdfplot([input.tDecisionTimeMs{:}]);
    set(pH, 'Color', 'g');
  if nLeft>0,
    pH1 = cdfplot([input.tDecisionTimeMs{leftTrialIx == 1}]);
      set(pH1, 'Color', 'b');
  end
  if nRight>0,
    pH2 = cdfplot([input.tDecisionTimeMs{leftTrialIx == 0}]);
      set(pH2, 'Color', [0.8 0.8 0]);
  set(gca, 'XLim', [0 decMax], ...
          'YLim', [0 1], ...
          'XTick', [0:decisionInterval:decMax]);
  end
elseif sum(block2Ix)>0
  pH1 = cdfplot([input.tDecisionTimeMs{block2Ix == 0 & leftTrialIx == 1}]);
    set(pH1, 'LineWidth', 2);
  hold on
  pH2 = cdfplot([input.tDecisionTimeMs{block2Ix == 0 & leftTrialIx == 0}]);
    set(pH2, 'LineStyle', '--', 'LineWidth', 2);
  pH3 = cdfplot([input.tDecisionTimeMs{block2Ix == 1 & leftTrialIx == 1}]);
    set(pH3, 'Color', yColor, 'LineWidth', 2);
  pH4 = cdfplot([input.tDecisionTimeMs{block2Ix == 1 & leftTrialIx == 0}]);
    set(pH4, 'Color', yColor, 'LineStyle', '--', 'LineWidth', 2);
    set(gca, 'XLim', [0 decMax], ...
          'YLim', [0 1], ...
          'XTick', [0:decisionInterval:decMax]);
end
grid on;
hold on;
title('Decision Time CDF: y=right,b=left');
xlabel('Time');


%%%%%%%%%%%%%%%%
%% 10 - Right/Left Contrast Decision

axH  = subplot(spSz{:},10);
hold on

if input.gratingContrastDiffSPO > 10
  contrastDifferenceRight = chop(cell2mat(input.rightGratingContrast) - cell2mat(input.leftGratingContrast),2);
elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
  contrastDifferenceRight = chop(cell2mat(input.rightGratingContrast) - cell2mat(input.leftGratingContrast),2);
elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
  contrastDifferenceRight =chop(cell2mat(input.rightGratingContrast) ./ cell2mat(input.leftGratingContrast),2);
end


plotTrsB1 = contrastDifferenceRight((correctIx|incorrectIx)&~block2Ix);
nLevelsB1 = unique(plotTrsB1);
percentContCellB1 = cell(1,length(nLevelsB1));
for kk=1:length(nLevelsB1)
    valB1 = nLevelsB1(kk);
    valIxB1 = contrastDifferenceRight==valB1;
    totalNTrialsValB1 = sum(contrastDifferenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix));
    if min(contrastDifferenceRight) < 0
      if valB1>=0,
          rightNTrialsValB1 = sum(contrastDifferenceRight((valIxB1&correctIx)&~block2Ix));
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
      elseif valB1<0,
          rightNTrialsValB1 = sum(contrastDifferenceRight((valIxB1&incorrectIx)&~block2Ix));
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
      end
    else
        if valB1>=1,
            rightNTrialsValB1 = sum(contrastDifferenceRight(valIxB1&correctIx&~block2Ix));
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
        elseif valB1<1,
            rightNTrialsValB1 = sum(contrastDifferenceRight(valIxB1&incorrectIx&~block2Ix));
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
        end
    end
end

if sum(block2Ix)>0
  plotTrsB2 = contrastDifferenceRight((correctIx|incorrectIx)&block2Ix);
  nLevelsB2 = unique(plotTrsB2);
  percentContCellB2 = cell(1,length(nLevelsB2));
  for kk=1:length(nLevelsB2)
    valB2 = nLevelsB2(kk);
    valIxB2 = contrastDifferenceRight==valB2;
    totalNTrialsValB2 = sum(contrastDifferenceRight(valIxB2&(correctIx|incorrectIx)&block2Ix));
    if min(contrastDifferenceRight) < 0
      if valB2>=0,
          rightNTrialsValB2 = sum(contrastDifferenceRight((valIxB2&correctIx)&block2Ix));
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      elseif valB2<0,
          rightNTrialsValB2 = sum(contrastDifferenceRight((valIxB2&incorrectIx)&block2Ix));
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      end
    else
        if valB2>=1,
            rightNTrialsValB2 = sum(contrastDifferenceRight(valIxB2&correctIx&block2Ix));
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        elseif valB2<1,
            rightNTrialsValB2 = sum(contrastDifferenceRight(valIxB2&incorrectIx&block2Ix));
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        end
    end
  end
end

if input.doNoGo
  didNoGoIx = celleqel2mat_padded(input.didNoGo);
  didNoGoIx
  plotTrsNoGo = contrastDifferenceRight(noGoIx|didNoGoIx);
  nLevelsNoGo = unique(plotTrsNoGo);
  for kk=1:length(nLevelsNoGo)
    valNoGo = nLevelsNoGo(kk);
    valIx = contrastDifferenceRight==valNoGo;
    totalNTrialsVal = sum(valIx);
    totalNTrialsValNoGo = sum(valIx & didNoGoIx);
    percentContCellNoGo{kk} = totalNTrialsValNoGo/totalNTrialsVal;
  end
end

if min(contrastDifferenceRight) < 0
    minX = min(contrastDifferenceRight);
    maxX = max(contrastDifferenceRight);
    xLimm = [minX maxX];
else
    minX = min(contrastDifferenceRight,[],2);
    maxX = max(contrastDifferenceRight,[],2);
    xLimm = [minX maxX];
        if ~(xLimm(1)==0),
            xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
        else
            xL1 = [0 ceil(log10(xLimm(2)))];
        end
    xTickL = 10.^(xL1(1):1:xL1(2));
    xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
    xTLabelL = cellstr(num2str(xTickL(:)));
end

pH1 = plot(nLevelsB1, cell2mat(percentContCellB1), 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
if sum(block2Ix)>= 1
  pH2 = plot(nLevelsB2, cell2mat(percentContCellB2), 'Color', yColor, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end
if input.doNoGo
  pH3 = plot(nLevelsNoGo, cell2mat(percentContCellNoGo), 'Color', 'r', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end
if input.gratingContrastDiffSPO <= 100
  vH = plot([1 1],[0 1]);
else
  vH = plot([0 0],[0 1]);
end
set(vH, 'Color', 'g');
set(gca, 'XLim', [minX maxX], ...
         'YLim', [0 1]);
if min(contrastDifferenceRight) < 0
      xlabel('Contrast Difference (R-L)')
      set(gca, 'XTick', [-1:0.5:1], ...
               'YTick', [0:0.25:1],...
               'XGrid', 'on');
else
    xlabel('Contrast Difference (R/L)')
    set(gca,'XScale', 'log', ...
            'XGrid', 'on',...
            'XTick', xTickL,...
            'XTickLabel', xTLabelL);
end

ylabel('% Right')

if input.doNoGo
  ylabel('% Right / % NoGo')
end

grid on

maxCell = length(nLevelsB1);
lowMidCell = maxCell/2;
highMidCell = lowMidCell + 1;
A = cell2mat(percentContCellB1(maxCell));
B = cell2mat(percentContCellB1(1));
M = (A + B)/2;
Selectivity = A - B;
rndSelectivity = roundn(Selectivity, -3);
Bias = M - 0.5;
rndBias = roundn(Bias, -2);

if mod(highMidCell,1) == 0
  C = cell2mat(percentContCellB1(highMidCell));
  D = cell2mat(percentContCellB1(lowMidCell));
  M2 = (C + D)/2;
  InnerBias = M2 - 0.5;
  rndInnerBias = roundn(InnerBias, -2);
else
  InnerBias = NaN;
  rndInnerBias = NaN;
end


title({['Selectivity = ' num2str(rndSelectivity)], ['Intrinsic/Challenged Bias = ' num2str(rndBias) '/' num2str(rndInnerBias)]});
%%%%%%%%%%%%%%%%%

%% 11 - Target Contrast x %Correct

axH = subplot(spSz{:},11);
hold on;

if nCorr>0 && input.doTestRobot==0,
    contTargetV = cell2mat_padded(input.tGratingContrast).*100;
    plotTrsB1 = contTargetV((correctIx|incorrectIx)&~block2Ix);
    uqTargetB1 = unique(plotTrsB1);
    nLevelsB1 = length(uqTargetB1);
    percentCellB1 = cell(1,nLevelsB1);
    for jj=1:nLevelsB1
        valB1 = uqTargetB1(jj);
        valIxB1 = contTargetV==valB1;
        totalNTrialsValB1 = sum(contTargetV(valIxB1&(correctIx|incorrectIx)&~block2Ix));
        corrNTrialsValB1 = sum(contTargetV(valIxB1&correctIx&~block2Ix));
        percentCellB1{jj} = corrNTrialsValB1/totalNTrialsValB1;
    end
    if sum(block2Ix)>0
      plotTrsB2 = contTargetV((correctIx|incorrectIx)&block2Ix);
    uqTargetB2 = unique(plotTrsB2);
    nLevelsB2 = length(uqTargetB2);
    percentCellB2 = cell(1,nLevelsB2);
    for jj=1:nLevelsB2
        valB2 = uqTargetB2(jj);
        valIxB2 = contTargetV==valB2;
        totalNTrialsValB2 = sum(contTargetV(valIxB2&(correctIx|incorrectIx)&block2Ix));
        corrNTrialsValB2 = sum(contTargetV(valIxB2&correctIx&block2Ix));
        percentCellB2{jj} = corrNTrialsValB2/totalNTrialsValB2;
    end
    end
    pH1 = plot(uqTargetB1, cell2mat(percentCellB1));
    if sum(block2Ix)>0
      ph2 = plot(uqTargetB2, cell2mat(percentCellB2),'Color', yColor);
    end
    minX = 0;
    maxX = 100;
    xLimm = [minX maxX];
    if ~(xLimm(1)==0),
      xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
    else
      xL1 = [0 ceil(log10(xLimm(2)))];
    end
    xTickL = 10.^(xL1(1):1:xL1(2));
    xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
    xTLabelL = cellstr(num2str(xTickL(:)));
    set(pH, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    set(gca, 'YLim', [0 1], ...
           'XLim', [minX maxX], ...
           'XScale', 'log', ...
           'XGrid', 'on');
       set(gca, 'XTick', xTickL);
       set(gca, 'XTickLabel', xTLabelL);
       ylabel('Correct (%)')
       xlabel('Target Contrast')
       title('Percent Correct by Target Contrast')
end
%%%%%%%%%%%%%%%%
if input.doBlocks==1

%% Block Performance
%% Indexing and pre-analysis
tLeftTrial = leftTrialIx;
leftBlockN = cell2mat_padded(input.tNBlockLeftTrsCompleted);
rightBlockN = cell2mat_padded(input.tNBlockRightTrsCompleted);
as = struct;

%% Left Analysis
if nLeft>0,
  axH  = subplot(spSz{:},11);
LBlockValues = unique(leftBlockN);
for i = 1:length(LBlockValues),
    trIx = (leftBlockN==LBlockValues(i)) & tLeftTrial;
    nTrs = sum(trIx);
    as.LnTrs(i) = nTrs;
    as.leftSuccess(i) = sum(trIx & correctIx)/nTrs; 
    as.leftIncorrs(i) = sum(trIx & incorrectIx)/nTrs; 
    as.leftIgnore(i) = sum(trIx & ignoreIx)/nTrs; 
end
hold on
plot(LBlockValues, as.leftSuccess, 'k', LBlockValues, as.leftIncorrs, 'g', LBlockValues, as.leftIgnore, 'm')
xlabel('Block Position')
ylabel('Fraction of Outcomes')
title('Left Block Analysis')
xlim([1 length(LBlockValues)])
end



%% Right Analysis
if nRight>0,
axH  = subplot(spSz{:},12);
RBlockValues = unique(rightBlockN);
for i = 1:length(RBlockValues),
    trIx = (rightBlockN==RBlockValues(i)) & ~tLeftTrial;
    nTrs = sum(trIx);
    as.RnTrs(i) = nTrs;
    as.rightSuccess(i) = sum(trIx & correctIx)/nTrs; 
    as.rightIncorrs(i) = sum(trIx & incorrectIx)/nTrs; 
    as.rightIgnore(i) = sum(trIx & ignoreIx)/nTrs; 
end
hold on
plot(RBlockValues, as.rightSuccess, 'k', RBlockValues, as.rightIncorrs, 'g', RBlockValues, as.rightIgnore, 'm')
xlabel('Block Position')
ylabel('Fraction of Outcomes')
title('Right Block Analysis')
xlim([1 length(RBlockValues)])
end
end
%%%%%%%%%%%%%%%%

end
end

%%%%%%%%%%%%%%%%
%set(gcf, 'Visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfunctions

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfig_print(epParamsIn{:});
end

function saveFile(hObject, evendata, input)
    save(input.savedDataName, 'input')
end

function outM = subCell2PadVect(cellV, padChar, nMaxTrials)
if nargin < 2, padChar = NaN; end
if nargin < 3, nMaxTrials = length(cellV); end
isEIx = cellfun(@isempty, cellV);
outM = repmat(padChar, [1 nMaxTrials]);
outM(~isEIx) = cell2mat(cellV(~isEIx));
end
