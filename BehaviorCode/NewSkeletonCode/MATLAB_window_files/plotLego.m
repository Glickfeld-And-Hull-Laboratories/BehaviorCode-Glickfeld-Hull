function input = plotLego(data_struct, input)

% essential consts
name = 'Lego';
cs = lego_constants;
spSz = {3,3};

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
 case 'MaunsellMouse1'
  figPos = [1111 208 806 923];
 otherwise
  figPos = [930 140 780 905];
end
set(figH, 'Position', figPos);

%% set up arrays
if ~isfield(input, 'changedStr')
  input.changedStr = {};
end

% some data processing
nTrial = length(input.trialOutcomeCell);

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
leftTrialIx = cell2mat(input.tLeftTrial);

% Left bias indexing
left_correct_ind = find((double(leftTrialIx)+double(correctIx))==2);
left_incorrect_ind = intersect(find(leftTrialIx==0), find(incorrectIx==1));
tLeftResponse = zeros(size(correctIx));
tLeftResponse([left_correct_ind left_incorrect_ind]) = 1;
tLeftResponse(find(ignoreIx)) = NaN;
input.tLeftResponse = tLeftResponse;

%  TODO:  this should be automated in a loop, or better should not have to convert from cell at all
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;

% stimulus str            
stimStr = strcat(mat2str(input.gratingHeightDeg), ' x ', mat2str(input.gratingWidthDeg), ' deg, ', ...
    num2str(input.gratingSpatialFreqCPD), ' cpd, ', mat2str(input.gratingEccentricityDeg), ' deg');

if input.gratingSpeedDPS > 0,
    stimStr = strcat(stimStr, ', ', mat2str(input.gratingSpeedDPS), ' dps, ', ...
        mat2str(input.gratingStartingPhaseDeg), ' deg');
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
         
	t2H(1) = text(0.00, 0.8, {'Trials:', 'Correct:', 'Incorrect:', 'Missed:'});
	t2H(2) = text(0.35, 0.8, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
				sprintf('%d', nInc), sprintf('%d', nIg)});
	t2H(3) = text(0.54, 0.8, {' ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
				sprintf('%.0f%%', nInc / numTrials * 100.0), ...
				sprintf('%.0f%%', nIg / numTrials * ...
                                        100.0)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');


        tStr = sprintf(['Decision median:\n', ...
                        '   %5.1f ms\n'], ...
                       median(decisionV(correctIx)));
               
        text(0.8, 0.8, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left');
         
            
        tStr = sprintf( ['Reaction Time: \t%5.2f s;   ITI %d ms \n', ...
                         'Reward : %d ms \n', ...
                         'Timeouts (ign,inc):\t%4.1f, %4.1f s\n', ...
                         'trPer80: %s\n', ...
                         'Stim: %s \n'], ...
                        input.reactionTimeMs/1000, ...
                        input.itiTimeMs, ...
                        input.rewardTimeUs/1000, ...
                        input.ignoreTimeoutMs/1000, ...
                        input.incorrectTimeoutMs/1000, ...
                        trPer80Str, ...
                        stimStr);

        text(0.0, 0.35, tStr, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

        set(gcf, 'Visible', 'on');
        
%%%%%%%%%%%%%%%

%% 2 - smoothed perf curve
axH = subplot(spSz{:},2:3);
hold on;
lH2 = plot(smooth(double(correctIx), 100, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(smooth(double(ignoreIx), 100, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2);

lH4 = plot(smooth(double(incorrectIx), 100, smoothType));
  set(lH4, 'Color', 'g', ...
           'LineWidth', 2);
  anystack(lH4, 'bottom');

  anystack(lH3, 'bottom');


title('Outcome Plot');
ylabel('Pct Occurrence');
set(gca, 'YLim', [0 1]);

xlabel('Trials');
trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);

%%%%%%%%%%%%%%%%
  
%% List changed text
axH = subplot(spSz{:}, 4);
hold on;
set(axH, 'Visible', 'off');
    set(axH, 'OuterPosition', [0.02 0.49, 0.25, 0.2])
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
    
title('Bias Plot');
ylabel('Pct Left');
set(gca, 'YLim', [0 1]);

xlabel('Trials');
trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
%%%%%%%%%%%%%%%%

%% 7 - Constrast Difference x Correct reaction times plot
axH = subplot(spSz{:},7);
hold on;

contDiffV = cell2mat_padded(input.tGratingContrast) - cell2mat_padded(input.dGratingContrast);
corrDiffV = contDiffV(correctIx);
uqDiff = unique(corrDiffV)
nLevels = length(uqDiff);
corrDiffCell = cell(1,nLevels);
decTimes = cell2mat_padded(input.tDecisionTimeMs);
corrDecTimes = decTimes(correctIx);

for ii = 1:length(uqDiff),
    val = uqDiff(ii);
    tix = corrDiffV==val;
    corrDiffCell{ii} = nanmean(corrDecTimes(tix));
end
minD = min(cell2mat_padded(corrDiffCell));
maxD = max(cell2mat_padded(corrDiffCell));
minX = min(uqDiff);
maxX = max(uqDiff);
xLimm = [minX maxX];

pH = plot(uqDiff, cell2mat(corrDiffCell));
xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
xTickL = 10.^(xL1(1):1:xL1(2));
xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
xTLabelL = cellstr(num2str(xTickL(:)));
set(pH, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 9);
set(gca, 'YLim', [minD-100 maxD+100], ...
           'XLim', [0 1], ...
           'XScale', 'log', ...
           'XGrid', 'on');
set(gca, 'XTick', xTickL);
set(gca, 'XTickLabel', xTLabelL);
ylabel('Decision Time (ms)')
xlabel('Contrast Difference')

%% Add a save button
if ~exist(cs.behavPdfPath, 'file')
  mkdir(cs.behavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  cs.behavPdfPath, ...
                  datestr(now, 'yymmdd'), input.subjectNum);
epParams = { figNum, outName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
bH = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF figure : %s', outName), ...
               'Units', 'pixels', ...
               'Position', [5 5 650 20], ...
               'Callback', { @saveButtonCb, epParams });

%%%%%%%%%%%%%%%%
%set(gcf, 'Visible', 'on');

end			
