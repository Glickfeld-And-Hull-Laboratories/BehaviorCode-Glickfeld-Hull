function input = plotLego(data_struct, input)

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
block2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);

block2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
block1Ix = block2Ix==0;

if isfield(input,'isNoGo')
noGoIx = cell2mat_padded(input.isNoGo);
else
noGoIx = zeros(size(block2Ix));
end
zeroConIx = zeros(size(block2Ix));
if input.doZeroConTrials
  totCon = sum(celleqel2mat_padded(input.tGratingContrast)+celleqel2mat_padded(input.dGratingContrast),1);
  zeroConIx(find(totCon==0)) = 1;
end

if isfield(input,'doMask')
    maskIx = celleqel2mat_padded(input.tDoMask);
    type1MaskIx = maskIx;
    if isfield(input,'doType2Mask')
      type2MaskIx = celleqel2mat_padded(input.tDoType2Mask);
      type1MaskIx(find(type2MaskIx)) = 0;
    else
      type2MaskIx = zeros(size(block2Ix));
    end
else
    maskIx = zeros(size(block2Ix));
    type2MaskIx = zeros(size(block2Ix));
    type1MaskIx = zeros(size(block2Ix));
end

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
if isfield(input, 'leftTrPer80Level1')
  leftTrPer80Str = ['[' num2str([input.leftTrPer80Level1 input.leftTrPer80Level2 input.leftTrPer80Level3 input.leftTrPer80Level4 input.leftTrPer80Level5 input.leftTrPer80Level6 input.leftTrPer80Level7 input.leftTrPer80Level8]) ']' ];
else
  leftTrPer80Str = [];
end

leftTrialIx = celleqel2mat_padded(input.tLeftTrial);
leftTrialIx(find(noGoIx)) = 0;
leftTrialIx(find(zeroConIx)) = 0;
leftTrNs = find(leftTrialIx);
rightTrialIx = ~leftTrialIx;
rightTrialIx(find(noGoIx)) = 0;
rightTrialIx(find(zeroConIx)) = 0;
rightTrNs = find(rightTrialIx);
leftTrialIx = logical(leftTrialIx);
nLeft = sum(leftTrialIx);
nRight = sum(rightTrialIx);
noGoIx = logical(noGoIx);
nNoGo = sum(noGoIx);

leftOutcomes = input.trialOutcomeCell(leftTrialIx);
leftCorr = strcmp(leftOutcomes, 'success');
nLeftCorr = sum(leftCorr);
leftIgn = strcmp(leftOutcomes, 'ignore');
nLeftIgn = sum(leftIgn);
leftInc = strcmp(leftOutcomes, 'incorrect');
nLeftInc = sum(leftInc);

rightOutcomes = input.trialOutcomeCell(rightTrialIx);
rightCorr = strcmp(rightOutcomes, 'success');
nRightCorr = sum(rightCorr);
rightIgn = strcmp(rightOutcomes, 'ignore');
nRightIgn = sum(rightIgn);
rightInc = strcmp(rightOutcomes, 'incorrect');
nRightInc = sum(rightInc);

noGoOutcomes = input.trialOutcomeCell(noGoIx);
noGoCorr = strcmp(noGoOutcomes, 'success');
nNoGoCorr = sum(noGoCorr);
noGoInc = strcmp(noGoOutcomes, 'incorrect');
nNoGoInc = sum(noGoInc);
noGoCorrIx = NaN(size(noGoIx));
noGoInd = find(noGoIx);
noGoCorrIx(noGoInd(noGoCorr)) = 1;
noGoCorrIx(noGoInd(noGoInc)) = 0;

%Left bias indexing
left_correct_ind = find((double(leftTrialIx)+double(correctIx))==2);
tLeftResponse = celleqel2mat_padded(input.tLeftResponse);
if input.doOriDiscrim
  tLeftResponse(intersect(find(celleqel2mat_padded(input.tLeftTrial)),find(strcmp(input.trialOutcomeCell, 'success')))) = 1;
  tLeftResponse(intersect(find(celleqel2mat_padded(input.tLeftTrial)),find(strcmp(input.trialOutcomeCell, 'incorrect')))) = 0;
  tLeftResponse(intersect(find(celleqel2mat_padded(input.tLeftTrial)==0),find(strcmp(input.trialOutcomeCell, 'success')))) = 0;
  tLeftResponse(intersect(find(celleqel2mat_padded(input.tLeftTrial)==0),find(strcmp(input.trialOutcomeCell, 'incorrect')))) = 1;
end

tRightResponse = celleqel2mat_padded(input.tRightResponse);
tLeftNoGo = nan(size(tLeftResponse));
leftNoGo = intersect(find(noGoIx),find(tLeftResponse));
rightNoGo = intersect(find(noGoIx),find(tRightResponse));
tLeftNoGo(leftNoGo) = 1;
tLeftNoGo(rightNoGo) = 0;
tLeftResponse(find(ignoreIx)) = NaN;
tLeftResponse(find(noGoIx)) = NaN;

%Right decision time indexing
right_correct_ind = find((double(rightTrialIx)+double(correctIx))==2);

%  TODO:  this should be automated in a loop, or better should not have to convert from cell at all
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;

% stimulus str 
if isfield(input,'gratingMaxDiameterDeg')           
  stimStr = strcat(mat2str(input.gratingMaxDiameterDeg), ' deg, ', ...
    num2str(input.gratingSpatialFreqCPD), ' cpd, ', mat2str(input.gratingEccentricityDeg), ' deg');
else
  stimStr = strcat(mat2str(input.gratingWidthDeg), ' deg, ', ...
    num2str(input.gratingSpatialFreqCPD), ' cpd, ', mat2str(input.gratingEccentricityDeg), ' deg');
end

if input.gratingSpeedDPS > 0,
    stimStr = strcat(stimStr, ', ', mat2str(input.gratingSpeedDPS), ' dps, ', ...
        mat2str(input.gratingStartingPhaseDeg), ' deg');
end

% leftStr formatting
if isfield(input,'doLeftSeparateOdds')
  if input.doLeftSeparateOdds == 0
    leftStr = sprintf('leftTrPer80: Matched');
  else
    leftStr= sprintf('leftTrPer80: %s \n', leftTrPer80Str);
  end
else
  leftStr = '';
end
% blockStr formatting
if input.doBlocks==0,
  blockStr = sprintf('Blocks: off, probLeft: %1.2f, probSwitch: %1.2f \n', double(input.tStimProbAvgLeft{nTrials}), double(input.stimProbAvgSwitch));

elseif input.doBlocks==1,
  blockStr = sprintf('Blocks: on, L: %2.0f, R:%2.0f \n', double(input.blockLeftTrs), double(input.blockRightTrs));
end

if input.doBlock2==1,
  block2Str = 'Block2: on \n';
  b2TrPer80Str = ['[' num2str([input.block2TrPer80Level1 input.block2TrPer80Level2 input.block2TrPer80Level3 input.block2TrPer80Level4 input.block2TrPer80Level5 input.block2TrPer80Level6 input.block2TrPer80Level7 input.block2TrPer80Level8]) ']' ];
  if input.block2DoTrialLaser==1,
    block2Str = sprintf('Block2: on, Trial Laser Power:%5.2fmW \n b2TrPer80: %s \n', double(input.block2TrialLaserPowerMw), b2TrPer80Str);
  end
elseif input.doBlock2==0,
  block2Str = 'Block2: off \n';
end

if isfield (input,'doAdapt')
  if input.doAdapt
    block2Str = ['Block1: Adapt to ' num2str(input.adaptDirectionDeg) ' \nBlock2: Adapt to' num2str(input.adaptDirectionDeg + input.adaptDirectionStepDeg) '\n'];
  end
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
    t2H(4) = text(0.9, 0.8, {'NoGo', sprintf('%d', nNoGo), sprintf('%d', nNoGoCorr), ...
                sprintf('%d', nNoGoInc), sprintf(''), ...
                sprintf('%5.0f', median(decisionV(noGoInc)))});
    t2H(5) = text(1.15, 0.8, {'Total', sprintf('%d', nTrial), sprintf('%d', nCorr), ...
                sprintf('%d', nInc), sprintf('%d', nIg), ...
                sprintf('%5.0f', median(decisionV(correctIx)))});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');


%         tStr = sprintf(['Decision median:\n', ...
%                         '   %5.1f ms\n'], ...
%                        median(decisionV(correctIx)));
%                
%         text(0.8, 0.8, tStr, ...
%              'VerticalAlignment', 'top', ...
%              'HorizontalAlignment', 'left');
         
if ~isfield(input,'doSizeDiscrim')
  input.doSizeDiscrim = 0;
end
        
if input.doContrastDiscrim || input.doContrastDetect || input.doSizeDiscrim
  decisionThreshold = (input.gratingEccentricityDeg-input.gratingTargetEccentricity)./ input.feedbackMotionSensitivity;
elseif input.doOriDiscrim
  decisionThreshold = input.gratingMaxDirectionDiff./ input.feedbackMotionSensitivity;
end

if ~isfield(input, 'stimOnTimeMs')
input.stimOnTimeMs = input.reactionTimeMs;
end

if isfield(input,'gratingMaxDiff')
  input.gratingMaxContrastDiff = input.gratingMaxDiff;
end


if isfield(input,'doContrastDiscrim')
  if input.doContrastDiscrim
    taskStr = sprintf(['MaxCon: %d ; SPO: %2.2f ; MaxDiff: %d ; SPO: %4.2f \n'] ,...
      input.gratingMaxContrast, input.gratingContrastSPO, input.gratingMaxContrastDiff, input.gratingContrastDiffSPO);
  elseif input.doSizeDiscrim
    taskStr = sprintf(['MaxSize: %d ; SPO: %2.2f ; MaxDiff: %d ; SPO: %4.2f \n'] ,...
      input.gratingMaxDiameterDeg, input.gratingDiameterSPO, input.gratingMaxDiameterDiff, input.gratingDiameterDiffSPO);
  elseif input.doOriDiscrim
    taskStr = sprintf(['MaxOriDiff: %d ; SPO: %2.2f \n'] ,...
      input.gratingMaxDirectionDiff, input.gratingDirectionDiffSPO);  
  end
else
  taskStr = sprintf(['MaxCon: %d ; SPO: %2.2f ; MaxDiff: %d ; SPO: %4.2f \n'] ,...
      input.gratingMaxContrast, input.gratingContrastSPO, input.gratingMaxContrastDiff, input.gratingContrastDiffSPO);
end
        if input.doFeedbackMotion
            FB_str = 'On';
        else
            FB_str = 'Off';
        end
        tStr = sprintf( ['Decision Time: \t%5.2f s;   ITI %d ms \n', ...
                         'Stim On Time: \t%5d ms \n', ...
                         'Thresholds (L,R):\t%4.0f pulses\n', ...
                         'Too-Fast Time: %d ms;  Reward: %s ms \n', ...
                         'Timeouts (ign,inc):\t%4.1f, %4.1f s\n', ...
                         'trPer80: %s\n', ...
                         taskStr, ...
                         leftStr, ...
                         'Stim: %s \n', ...
                         blockStr, ...
                         block2Str, ...
                         ['Feedback' FB_str ', Motion Sensitivity: \t%2.3f \n']], ...
                        input.reactionTimeMs/1000, ...
                        input.itiTimeMs, ...
                        input.stimOnTimeMs, ...
                        decisionThreshold, ...
                        input.tooFastTimeMs, ...
                        rewStr, ...
                        input.ignoreTimeoutMs/1000, ...
                        input.incorrectTimeoutMs/1000, ...
                        trPer80Str, ...
                        stimStr, ...
                        input.feedbackMotionSensitivity);

        text(0.0, 0.15, tStr, ...
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

lH5 = plot(smooth(double(noGoCorrIx), nTrials/5, smoothType));
set(lH5, 'Color', 'b', ...
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

  text(0,0.5,input.changedStr, ...
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

leftTrialIxDouble = double(leftTrialIx);
leftTrialIxDouble(find(noGoIx)) = NaN;
lh1 = plot(smooth(leftTrialIxDouble, amtSmooth, smoothType));
set(lh1, 'Color', 'k', 'LineWidth', 2);

lh2 = plot(smooth(double(tLeftResponse), amtSmooth, smoothType));
set(lh2,'Color', 'r', 'LineWidth', 2);

lh3 = plot(smooth(double(tLeftNoGo), amtSmooth, smoothType));
set(lh3,'Color', 'b', 'LineWidth', 2);

lh4= refline(0,0.5);
set(lh4, 'LineStyle', '--');
    
title('Bias Plot: Red = Resp, Black = Pres, Blue = NoGo');
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
decTimes = celleqel2mat_padded(input.tDecisionTimeMs);
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

%hack for switching between detect and discrim
dGratingContrast = celleqel2mat_padded(input.dGratingContrast);
ind = find(dGratingContrast == 0);
if length(ind)>0
  dGratingContrast(ind) = 0.0001;
end

if nCorr>0  %&& input.doTestRobot==0,
    if input.doContrastDiscrim
      if input.doContrastDetect == 0
        contDiffV = chop(celleqel2mat_padded(input.tGratingContrast) ./ dGratingContrast,2);
      else
        contDiffV = chop(celleqel2mat_padded(input.tGratingContrast) - dGratingContrast,2);
      end
    elseif input.doSizeDiscrim
      if input.gratingDiameterDiffSPO<10
        contDiffV = chop(celleqel2mat_padded(input.tGratingDiameterDeg) ./ celleqel2mat_padded(input.dGratingDiameterDeg),2);
      else
        contDiffV = chop(celleqel2mat_padded(input.tGratingDiameterDeg) - celleqel2mat_padded(input.dGratingDiameterDeg),2);
      end
    elseif input.doOriDiscrim
      contDiffV = chop(celleqel2mat_padded(input.tGratingDirectionStart) - double(input.gratingTargetDirection),2);
    end


    corrDiffV = contDiffV(correctIx);
    uqDiff = unique(corrDiffV);
    nLevels = length(uqDiff);
    corrDiffCell = cell(1,nLevels);
    decTimes = celleqel2mat_padded(input.tDecisionTimeMs);
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
    if isfield(input,'doContrastDiscrim')
      if input.doContrastDiscrim
        if isfield(input, 'dGratingContrastDiff')
          possDiffV = double(input.gratingMaxContrastDiff) ./ (2 .^ (lev./double(input.gratingContrastDiffSPO)))+1;

          minX = min(possDiffV,[],2);
          maxX = max(possDiffV,[],2);
        else
          possDiffV = input.gratingMaxContrastDiff ./ (2 .^ (lev./input.gratingContrastSPO));
          minX = min(possDiffV,[],2);
          maxX = max(possDiffV,[],2);
        end
      elseif input.doSizeDiscrim
        possDiffV = double(input.gratingMaxDiameterDiff) ./ (2 .^ (lev./double(input.gratingDiameterDiffSPO)))+1;
        minX = min(possDiffV,[],2);
        maxX = max(possDiffV,[],2);
      elseif input.doOriDiscrim
        possDiffV = double(input.gratingMaxDirectionDiff) ./ (2 .^ (lev./double(input.gratingDirectionDiffSPO)))+1;
        minX = min(possDiffV,[],2);
        maxX = max(possDiffV,[],2);
      end
    else
      if isfield(input, 'dGratingContrastDiff')
        possDiffV = double(input.gratingMaxContrastDiff) ./ (2 .^ (lev./double(input.gratingContrastDiffSPO)))+1;

        minX = min(possDiffV,[],2);
        maxX = max(possDiffV,[],2);
      else
        possDiffV = input.gratingMaxContrastDiff ./ (2 .^ (lev./input.gratingContrastSPO));
        minX = min(possDiffV,[],2);
        maxX = max(possDiffV,[],2);
      end
    end


    minD = min(cell2mat(corrDiffCell));
    maxD = max(cell2mat(corrDiffCell));
    
    hold on
    xLimm = [minX maxX];
    uqDiff = uqDiff;
    
    sc_num = 1; 
    
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
           'XGrid', 'on');


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
         if input.doContrastDiscrim
           if input.doContrastDetect == 0
             xlabel('Contrast Ratio (T/D)')
             title('Decision Time by Contrast Ratio')
           else
             xlabel('Target Contrast')
             title('Decision Time by Target Contrast')
            end
        elseif input.doSizeDiscrim
          if input.gratingDiameterDiffSPO<10
             xlabel('Size Difference (T/D)')
           else
             xlabel('Size Difference (T-D)')
            end
           title('Decision Time by Size Ratio')
        elseif input.doOriDiscrim
           xlabel('Ori Difference')
           title('Decision Time by Ori Difference')
        end
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
    didNoGoIx = celleqel2mat_padded(input.didNoGo);
    plotTrsNoGo = contDiffV(~noGoIx&didNoGoIx);
    if find(didNoGoIx)
      nLevelsNoGo = plotTrsNoGo;
      for kk=1:length(nLevelsNoGo)
        valNoGo = nLevelsNoGo(kk);
        valIx = contDiffV==valNoGo;
        totalNTrialsVal = sum(valIx);
        totalNTrialsValNoGo = sum(valIx & didNoGoIx);
        percentContCellNoGoByDiff{kk} = totalNTrialsValNoGo/totalNTrialsVal;
      end
    end
    if input.doOriDiscrim == 0 
    if ~find(uqDiff <=0)
      pH = plot(uqDiff, cell2mat(percentCell));
    else
      percentCell = cell2mat(percentCell);
      percentCell(find(uqDiff<=0)) = [];
      uqDiff(find(uqDiff<=0)) = [];
      pH = plot(uqDiff, percentCell);
    end
    if length(plotTrsNoGo)>0
       pH1 = plot(sc_num * nLevelsNoGo, cell2mat(percentContCellNoGoByDiff), 'Color', 'r');
       set(pH1, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    end

    set(pH, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    ylabel('Correct/Ignore (%)')
    set(gca, 'YLim', [0 1], ...
           'XGrid', 'on',...
           'Xscale', 'log');
       

       if isfield(input,'doContrastDiscrim')
         if input.doContrastDiscrim
           if input.doContrastDetect == 0
             xlabel('Contrast Ratio (T/D)')
             title('Percent Correct by Contrast Ratio')
           else
             xlabel('Target Contrast')
             title('Percent Correct by Target Contrast')
           end
        elseif input.doSizeDiscrim
          if input.gratingDiameterDiffSPO<10
            xlabel('Size Ratio (T/D)')
           else
             xlabel('Size Difference (T-D)')
            end
           title('Percent Correct by Size Ratio')
        elseif input.doOriDiscrim
           xlabel('Ori Difference')
           title('Percent Correct by Ori Difference')
      else
        if input.doContrastDetect == 0
         xlabel('Contrast Ratio (T/D)')
        else
         xlabel('Contrast Difference (T-D)')
        end
      title('Percent Correct by Contrast RAtio')
         end
    end
    end
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

if find(isnan(block2Ix))
  block2Ix(find(isnan(block2Ix))) = 0;
end

if sum(block2Ix)== 0
  pH = cdfplot([input.tDecisionTimeMs{~noGoIx}]);
    set(pH, 'Color', 'g');
  if nLeft>0,
    pH1 = cdfplot([input.tDecisionTimeMs{leftTrialIx}]);
      set(pH1, 'Color', 'b');
  end
  if nRight>0,
    pH2 = cdfplot([input.tDecisionTimeMs{rightTrialIx}]);
      set(pH2, 'Color', [0.8 0.8 0]);
  set(gca, 'XLim', [0 decMax], ...
          'YLim', [0 1], ...
          'XTick', [0:decisionInterval:decMax]);
  end
  if nNoGo>0,
    pH1 = cdfplot([input.tDecisionTimeMs{noGoIx}]);
      set(pH1, 'Color', 'c');
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

  if input.doContrastDiscrim
    if input.doContrastDetect
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif input.gratingContrastDiffSPO > 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
      differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
    end
  elseif input.doSizeDiscrim
    if input.gratingDiameterDiffSPO > 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingDiameterDeg) - celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    elseif ~isfield(input, 'dGratingDiameterDiff') & input.gratingDiameterDiffSPO <= 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingDiameterDeg) - celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    elseif isfield(input, 'dGratingDiameterDiff') & input.gratingDiameterDiffSPO <= 10
        differenceRight =chop(celleqel2mat_padded(input.rightGratingDiameterDeg) ./ celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    end
  elseif input.doOriDiscrim
      differenceRight = chop(celleqel2mat_padded(input.tGratingDirectionStart) - double(input.gratingTargetDirection),2);
      if isfield(input,'doMask')
          differenceRight_mask = chop(celleqel2mat_padded(input.tPlaidDirectionStart) - double(input.gratingTargetDirection),2);
      else
          differenceRight_mask = nan(size(differenceRight));
      end
  else
  if input.gratingContrastDiffSPO > 10
    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
  elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
  elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
    differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
  end
end

plotTrsB1 = differenceRight((correctIx|incorrectIx)&~block2Ix&~maskIx);
nLevelsB1 = unique(plotTrsB1);
curr100 = zeros(size(block2Ix));
nT= size(curr100,2);
if nT>100
  curr100(nT-99:nT) = 1;
end
percentContCellB1 = cell(1,length(nLevelsB1));
percentContCellB1_100 = cell(1,length(nLevelsB1));
for kk=1:length(nLevelsB1)
    valB1 = nLevelsB1(kk);
    valIxB1 = differenceRight==valB1;
    totalNTrialsValB1 = length(differenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix&~maskIx));
    totalNTrialsValB1_100 = length(differenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix&curr100&~maskIx));
    if min(differenceRight) < 0
      if valB1>=0,
          ind = setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]);  
          rightNTrialsValB1 = length(ind);
          ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]));  
          rightNTrialsValB1_100 = length(ind2);
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
          percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
      elseif valB1<0,
          ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]);
          rightNTrialsValB1 = length(ind);
          ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]));
          rightNTrialsValB1_100 = length(ind2);
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
          percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
      end
    else
        if valB1>=1,
            ind = setdiff(intersect(find(valIxB1),find(correctIx)), find(block2Ix));
            rightNTrialsValB1 = length(ind);
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
            ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]));
            rightNTrialsValB1_100 = length(ind2);
            percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
        elseif valB1<1,
            ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), find(block2Ix));
            rightNTrialsValB1 = length(ind);
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
            ind2 = intersect(find(curr100),setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]));
            rightNTrialsValB1_100 = length(ind2);
            percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
        end
    end
end

if sum(maskIx)>0
    if sum(type1MaskIx)>0
      plotTrsM1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type1MaskIx);
      nLevelsM1 = unique(plotTrsM1);
      percentContCellM1 = cell(1,length(nLevelsM1));
      for kk=1:length(nLevelsM1)
          valM1 = nLevelsM1(kk);
          valIxM1 = differenceRight_mask==valM1;
          totalNTrialsValM1 = length(differenceRight_mask(valIxM1&(correctIx|incorrectIx)&~block2Ix&type1MaskIx));
          if min(differenceRight_mask) < 0
            if valM1>=0,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                rightNTrialsValM1 = length(ind);
                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
            elseif valM1<0,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                rightNTrialsValM1 = length(ind);
                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
            end
          else
              if valM1>=1,
                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                  rightNTrialsValM1 = length(ind);
                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
              elseif valM1<1,
                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                  rightNTrialsValM1 = length(ind);
                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
              end
          end
      end
    end
    if sum(type2MaskIx)>0
      plotTrsT2M1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type2MaskIx);
      nLevelsT2M1 = unique(plotTrsT2M1);
      percentContCellT2M1 = cell(1,length(nLevelsT2M1));
      for kk=1:length(nLevelsT2M1)
          valT2M1 = nLevelsT2M1(kk);
          valIxT2M1 = differenceRight_mask==valT2M1;
          totalNTrialsValT2M1 = length(differenceRight_mask(valIxT2M1&(correctIx|incorrectIx)&~block2Ix&type2MaskIx));
          if min(differenceRight_mask) < 0
            if valT2M1>=0,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                rightNTrialsValT2M1 = length(ind);
                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
            elseif valT2M1<0,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                rightNTrialsValT2M1 = length(ind);
                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
            end
          else
              if valT2M1>=1,
                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                  rightNTrialsValT2M1 = length(ind);
                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
              elseif valT2M1<1,
                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                  rightNTrialsValT2M1 = length(ind);
                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
              end
          end
      end
    end
end

if sum(block2Ix)>0
  plotTrsB2 = differenceRight((correctIx|incorrectIx)&block2Ix);
  nLevelsB2 = unique(plotTrsB2);
  percentContCellB2 = cell(1,length(nLevelsB2));
  for kk=1:length(nLevelsB2)
    valB2 = nLevelsB2(kk);
    valIxB2 = differenceRight==valB2;
    totalNTrialsValB2 = length(differenceRight(valIxB2&(correctIx|incorrectIx)&block2Ix));
    if min(differenceRight) < 0
      if valB2>=0,
          ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
          rightNTrialsValB2 = length(ind);
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      elseif valB2<0,
          ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
          rightNTrialsValB2 = length(ind);
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      end
    else
        if valB2>=1,
            ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
            rightNTrialsValB2 = length(ind);
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        elseif valB2<1,
            ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
            rightNTrialsValB2 = length(ind);
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        end
    end
  end
end

if sum(maskIx&block2Ix)>0
  if sum(type1MaskIx)>0 type1MaskIx
    plotTrsM2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type1MaskIx);
    nLevelsM2 = unique(plotTrsM2);
    percentContCellM2 = cell(1,length(nLevelsM2));
    for kk=1:length(nLevelsM2)
        valM2 = nLevelsM2(kk);
        valIxM2 = differenceRight_mask==valM2;
        totalNTrialsValM2 = length(differenceRight_mask(valIxM2&(correctIx|incorrectIx)&block2Ix&type1MaskIx));
        if min(differenceRight_mask) < 0
          if valM2>=0,
              ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
              rightNTrialsValM2 = length(ind);
              percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
          elseif valM2<0,
              ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
              rightNTrialsValM2 = length(ind);
              percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
          end
        else
            if valM2>=1,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
                rightNTrialsValM2 = length(ind);
                percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
            elseif valM2<1,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
                rightNTrialsValM2 = length(ind);
                percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
            end
        end
    end
  end
  if sum(type2MaskIx)>0
    plotTrsT2M2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type2MaskIx);
    nLevelsT2M2 = unique(plotTrsT2M2);
    percentContCellT2M2 = cell(1,length(nLevelsT2M2));
    for kk=1:length(nLevelsT2M2)
        valT2M2 = nLevelsT2M2(kk);
        valIxT2M2 = differenceRight_mask==valT2M2;
        totalNTrialsValT2M2 = length(differenceRight_mask(valIxT2M2&(correctIx|incorrectIx)&block2Ix&type2MaskIx));
        if min(differenceRight_mask) < 0
          if valT2M2>=0,
              ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
              rightNTrialsValT2M2 = length(ind);
              percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
          elseif valT2M2<0,
              ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
              rightNTrialsValT2M2 = length(ind);
              percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
          end
        else
            if valT2M2>=1,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
                rightNTrialsValT2M2 = length(ind);
                percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
            elseif valT2M2<1,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
                rightNTrialsValT2M2 = length(ind);
                percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
            end
        end
    end
  end
end


didNoGoIx = celleqel2mat_padded(input.didNoGo);
plotTrsNoGo = differenceRight(noGoIx|didNoGoIx);
nLevelsNoGo = unique(plotTrsNoGo);
for kk=1:length(nLevelsNoGo)
  valNoGo = nLevelsNoGo(kk);
  valIx = differenceRight==valNoGo;
  totalNTrialsVal = sum(valIx);
  totalNTrialsValNoGo = sum(valIx & didNoGoIx);
  percentContCellNoGo{kk} = totalNTrialsValNoGo/totalNTrialsVal;
end

if min(differenceRight) < 0
    minX = min(differenceRight);
    maxX = max(differenceRight);
    xLimm = [minX maxX];
else
    minX = min(differenceRight,[],2);
    maxX = max(differenceRight,[],2);
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

nLevelsB1 = nLevelsB1(~isnan(nLevelsB1));
pH1 = plot(nLevelsB1, cell2mat(percentContCellB1), 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
pH1x = plot(nLevelsB1, cell2mat(percentContCellB1_100), 'k', 'LineWidth', .5);
if sum(block2Ix)>= 1
  pH2 = plot(nLevelsB2, cell2mat(percentContCellB2), 'Color', yColor, 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end

if sum(didNoGoIx)>0
  pH3 = plot(nLevelsNoGo, cell2mat(percentContCellNoGo), 'Color', 'r', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end

if sum(type1MaskIx)>= 1
  pH4 = plot(nLevelsM1, cell2mat(percentContCellM1), 'Color', 'c', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end
if sum(type1MaskIx&block2Ix)>= 1
  pH5 = plot(nLevelsM2, cell2mat(percentContCellM2), 'Color', '--c', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end

if sum(type2MaskIx)>= 1
  pH6 = plot(nLevelsT2M1, cell2mat(percentContCellT2M1), 'Color', 'g', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end
if sum(type2MaskIx&block2Ix)>= 1
  pH7 = plot(nLevelsT2M2, cell2mat(percentContCellT2M2), 'Color', '--g', 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
end


if input.gratingContrastDiffSPO <= 100
  vH = plot([1 1],[0 1]);
else
  vH = plot([0 0],[0 1]);
end
set(vH, 'Color', 'g');

set(gca, 'XLim', [minX maxX], ...
         'YLim', [0 1]);
if min(differenceRight) < 0
  if isfield(input,'doContrastDiscrim')
    if input.doContrastDiscrim
        xlabel('Contrast Difference (R-L)')
        set(gca, 'XTick', [-1:0.5:1], ...
                 'YTick', [0:0.25:1],...
                 'XGrid', 'on');
    elseif input.doSizeDiscrim
        xlabel('Size Difference (R-L)')
        set(gca, 'XTick', [-input.gratingMaxDiameterDeg:5:input.gratingMaxDiameterDeg], ...
                 'YTick', [0:0.25:1],...
                 'XGrid', 'on');
    elseif input.doOriDiscrim
        xlabel('Ori Difference')
        set(gca, 'XTick', [-input.gratingMaxDirectionDiff:5:input.gratingMaxDirectionDiff], ...
                 'YTick', [0:0.25:1],...
                 'XGrid', 'on');
    end
  else
    xlabel('Contrast Difference (R-L)')
        set(gca, 'XTick', [-1:0.5:1], ...
                 'YTick', [0:0.25:1],...
                 'XGrid', 'on');
  end
else
  set(gca,'XScale', 'log', ...
            'XGrid', 'on',...
            'XTick', xTickL,...
            'XTickLabel', xTLabelL);
  if isfield(input,'doContrastDiscrim')
    if input.doContrastDiscrim
      xlabel('Contrast Difference (R/L)')
    elseif input.doSizeDiscrim
      xlabel('Size Difference (R/L)')
    elseif input.doOriDiscrim
      xlabel('Ori Difference')
    end
  else
    xlabel('Contrast Difference (R/L)')
  end
end

if input.doOriDiscrim
    ylabel('% Left/Ignore')
else
    ylabel('% Right/Ignore')
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

if nCorr>0 %&& input.doTestRobot==0,
  if isfield(input,'doContrastDiscrim')
    if input.doContrastDiscrim
      contTargetV = celleqel2mat_padded(input.tGratingContrast);
    elseif input.doSizeDiscrim
      contTargetV = celleqel2mat_padded(input.tGratingDiameterDeg);
    elseif input.doOriDiscrim
      contTargetV = abs(celleqel2mat_padded(input.tGratingDirectionStart)-double(input.gratingTargetDirection));
    end
  else
    contTargetV = celleqel2mat_padded(input.tGratingContrast).*100;
  end
    plotTrsB1 = contTargetV((correctIx|incorrectIx)&~block2Ix);
    uqTargetB1 = unique(plotTrsB1);
    nLevelsB1 = length(uqTargetB1);
    percentCellB1 = cell(1,nLevelsB1);
    for jj=1:nLevelsB1
        valB1 = uqTargetB1(jj);
        valIxB1 = contTargetV==valB1;
        totalNTrialsValB1 = length(contTargetV(valIxB1&(correctIx|incorrectIx)&~block2Ix));
        corrNTrialsValB1 = length(contTargetV(valIxB1&correctIx&~block2Ix));
        percentCellB1{jj} = corrNTrialsValB1/totalNTrialsValB1;
    end

    didNoGoIx = celleqel2mat_padded(input.didNoGo);
    plotTrsNoGo = contTargetV(~noGoIx&didNoGoIx);
    nLevelsNoGo = unique(plotTrsNoGo);
    for kk=1:length(nLevelsNoGo)
      valNoGo = nLevelsNoGo(kk);
      valIx = contTargetV==valNoGo;
      totalNTrialsVal = sum(valIx);
      totalNTrialsValNoGo = sum(valIx & didNoGoIx);
      percentContCellNoGoByTarg{kk} = totalNTrialsValNoGo/totalNTrialsVal;
    end
    if sum(block2Ix)>0
      plotTrsB2 = contTargetV((correctIx|incorrectIx)&block2Ix);
      uqTargetB2 = unique(plotTrsB2);
      nLevelsB2 = length(uqTargetB2);
      percentCellB2 = cell(1,nLevelsB2);
      for jj=1:nLevelsB2
          valB2 = uqTargetB2(jj);
          valIxB2 = contTargetV==valB2;
          totalNTrialsValB2 = length(contTargetV(valIxB2&(correctIx|incorrectIx)&block2Ix));
          corrNTrialsValB2 = length(contTargetV(valIxB2&correctIx&block2Ix));
          percentCellB2{jj} = corrNTrialsValB2/totalNTrialsValB2;
      end
    end

    pH1 = plot(uqTargetB1, cell2mat(percentCellB1));
    if sum(block2Ix)>0
      ph2 = plot(uqTargetB2, cell2mat(percentCellB2),'Color', yColor);
        set(ph2, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    end
    if length(plotTrsNoGo)>0
      ph3 = plot(nLevelsNoGo, cell2mat(percentContCellNoGoByTarg),'Color', 'r');
      set(ph3, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
    end
    if input.doContrastDiscrim
      minX = 0.01;
      maxX = 1;
    elseif input.doSizeDiscrim
      minX = 1;
      maxX = 100;
    elseif input.doOriDiscrim
      minX = 1;
      maxX = 45;
    end
    xLimm = [minX maxX];
    if ~(xLimm(1)==0),
      xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
    else
      xL1 = [0 ceil(log10(xLimm(2)))];
    end
    xTickL = 10.^(xL1(1):1:xL1(2));
    xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
    xTLabelL = cellstr(num2str(xTickL(:)));
    set(pH1, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);

    set(gca, 'YLim', [0 1], ...
           'XLim', [minX maxX], ...
           'XScale', 'log', ...
           'XGrid', 'on');
       ylabel('Correct/Ignore (%)')
       if isfield(input,'doContrastDiscrim')
         if input.doContrastDiscrim
           xlabel('Target Contrast')
           title('Percent Correct by Target Contrast')
         elseif input.doSizeDiscrim
           xlabel('Target Size')
           title('Percent Correct by Target Size')
         elseif input.doOriDiscrim
           xlabel('Abs Ori diff')
           title('Percent Correct by Ori Diff')
         end
        else
          xlabel('Target Contrast')
          title('Percent Correct by Target Contrast')
        end
end
%%%%%%%%%%%%%%%%
if input.doBlocks==1

%% Block Performance
%% Indexing and pre-analysis
tLeftTrial = leftTrialIx;
leftBlockN = celleqel2mat_padded(input.tNBlockLeftTrsCompleted);
rightBlockN = celleqel2mat_padded(input.tNBlockRightTrsCompleted);
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

%% Add a save button
if ~exist(cs.behavPdfPath, 'file')
  mkdir(cs.behavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  cs.behavPdfPath, ...
                  datestr(now, 'yymmdd-HHMM'), input.subjectNum);
dbName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  cs.dbPath, ...
                  datestr(now, 'yymmdd-HHMM'), input.subjectNum);
              
matParams = {input.savedDataName, 'input'};
epParams = { figNum, outName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
dbParams = { figNum, dbName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
bH1 = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF Figure Locally'), ...
               'Units', 'pixels', ...
               'Position', [5 5 365 20], ...
               'Callback', { @saveButtonCb, epParams});

bH2 = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf('SAVE DATA FILE'), ...
               'Units', 'pixels', ...
              'Position', [465 5 350 20], ...
               'Callback', {@saveFile,input});
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
