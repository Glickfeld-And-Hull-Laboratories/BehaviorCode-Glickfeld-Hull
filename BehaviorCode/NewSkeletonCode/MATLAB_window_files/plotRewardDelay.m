function input = plotRewardDelay(data_struct, input)

% essential consts
figNum = 4;
name = 'RewardDelay';
consts = reward_delay_constants;
subplotSz = {3,3};

smoothType = 'lowess';
%% draw figure
figH = figure(figNum);
set(figH, 'ToolBar', 'none');
clf;
figPos = [930 140 905 905];
set(figH, 'Position', figPos);
orange = [1 0.64 0];

%% set up arrays
if ~isfield(input, 'changedStr')
  input.changedStr = {};
end

% some data processing
nPts = length(input.holdTimesMs);
nTrial = length(input.trialOutcomeCell);
reqHoldStartV = celleqel2mat_padded(input.tReqHoldToStartMs);
successIx = strcmp(input.trialOutcomeCell, 'success');
earlyIx = strcmp(input.trialOutcomeCell, 'early');
lateIx = strcmp(input.trialOutcomeCell, 'late');
noreleaseIx = strcmp(input.trialOutcomeCell, 'norelease');
nCorr = sum(successIx);
nFail = sum(earlyIx);
nIg = sum(lateIx);
nNR = sum(noreleaseIx);
holdStarts = [input.holdStartsMs{:}];
nTrials = length(input.trialOutcomeCell);

reactV = celleqel2mat_padded(input.reactTimesMs);
holdV = celleqel2mat_padded(input.holdTimesMs);  % sometimes on client restart have empty elements here; should figure out why
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;
tTrialN = input.trialSinceReset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Values

if length(holdStarts) > 2
	holdTimeMS = input.holdTimesMs{end};
	switch input.trialOutcomeCell{end}
	case 'success'
		outcomeString = 'correct';
	case 'early'
		outcomeString = 'early';
	case 'late'
		outcomeString = 'late';
    case 'norelease'
        outcomeString = 'norelease';  
	end

        axH = subplot(subplotSz{:}, 1);						% default axes are 0 to 1
        
	set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

        numTrials = nTrial;
        
        set(gcf, 'Visible', 'off'); % hide figure during text
                                    % drawing - kludge

        text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 18);

        elMin = round((now - datenum(input.startDateVec)) * 24*60);
        startStr = datestr(input.startDateVec, 'HH:MM');
        text(0.00, 1.0, {'Subject Number:', 'Start time + elapsed:', 'Reward vol (m\pm std):'}, 'FontSize', 12);
	text(0.70, 1.0, ...
             { sprintf('%2d', input.subjectNum), ...
               sprintf('%s + %2dm', ...
                       startStr, elMin), ...
               sprintf('%.1f s \t(%g ms\\pm %g ms)', ...
                       nansum(juiceTimesMsV./1000), ...
                       chop(nanmean(juiceTimesMsV),2), ...
                       chop(nanstd(juiceTimesMsV),2)), ...
             }, 'FontSize', 12);
	t2H(1) = text(0.00, 0.8, {'Trials:', 'Correct:', 'Early:', 'Late:', 'No Release:'});
	t2H(2) = text(0.38, 0.8, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
				sprintf('%d', nFail), sprintf('%d', nIg), sprintf('%d', nNR)});
	t2H(3) = text(0.54, 0.8, {' ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
				sprintf('%.0f%%', nFail / numTrials * 100.0), ...
				sprintf('%.0f%%', nIg / numTrials *100.0), ...
                sprintf('%.0f%%', nNR / numTrials *100.0)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left', 'FontSize', 12);


        tStr = sprintf(['Hold, react median:\n', ...
                        '   %5.1f ms, %5.1f ms\n', ...
                        '\n', ...
                        '%s'], ...
                       nanmedian(holdV), ...
                       nanmedian(reactV(successIx)));
               
        text(0.8, 0.8, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);

        
        tStr = sprintf( ['Hold Time To Start Trial: \t%3d ms \n', ...
                         'Delay Time: \t%3d ms \n', ...
                         'WaitForLate Time: \t%3d ms \n', ...
                         'Timeouts (e,l,nr):\t%4.1fs, %4.1fs, %4.1fs\n', ...
                         'Reward Window Width:\t%5.0f ms \n',...
                         'Inter-Trial Interval: %d ms\n', ...
                         'Reward min, max: %3.0f, %3.0f ms   %s\n'], ...
                        input.reqHoldToStartMs, ...
                        input.delayTimeMs, ...
                        input.waitForLateTimeMs, ...
                        input.earlyTimeoutMs/1000, ...
                        input.lateTimeoutMs/1000, ...
                        input.noReleaseTimeoutMs/1000, ...
                        input.rewardWindowWidthMs, ...
                        input.itiTimeMs, ...
                        input.minRewardUs/1000, ...
                        input.maxRewardUs/1000);

        text(0.0, 0.3, tStr, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top', 'FontSize', 12);

        set(gcf, 'Visible', 'on');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Hold time histogram

axH = subplot(subplotSz{:}, 5);
if input.delayTimeMs  < 500
  maxX = 2000;
  disp('true');
else
  maxFail = double(input.delayTimeMs+input.reqHoldToStartMs+input.waitForLateTimeMs);
  maxX = ceil((maxFail)./500)*500;  % round up to nearest 500 ms.
end
visIx = holdV <= maxX;
nVisPts = sum(visIx);
if nVisPts > 50
  binWidth = iqr(holdV(visIx)) ./ nVisPts.^(1/3);	% robust version of std
  nBins = ceil(maxX ./ binWidth);
else
  nBins = 10;
end
edges = linspace(0, maxX, nBins);
Ns = histc(holdV(find(successIx)), edges);
Nf = histc(holdV(find(earlyIx)), edges);
Nl = histc(holdV(find(lateIx)), edges);
Nnr = histc(holdV(find(noreleaseIx)), edges);
if sum(Ns) + sum(Nf) + sum(Nl) + sum(Nnr) > 0
  bH = bar(edges, [Ns(:) Nf(:) Nl(:) Nnr(:)], 'stacked');%, edges, Nf, 'c', edges, Nl, 'm');
  set(bH, 'BarWidth', 1, 'LineStyle', 'none');
  set(bH(1),'facecolor','g')
  set(bH(2),'facecolor','c')
  set(bH(3),'facecolor','m')
  set(bH(4), 'facecolor', orange)
end
hold on;
xLim = [0 maxX];
set(gca, 'XLim', xLim);
yLim = get(gca, 'YLim');
if (length(get(gca, 'XTick')) > 4)
  xT = (0:500:maxX);
  set(gca, 'XTick', xT);
end

if ~isempty(Nf)
  pH = plot(edges, Nf);
  set(pH, 'LineStyle', '--', ...
          'Color', 'c');
end
if ~isempty(Nl)
  pH = plot(edges, Nl);
  set(pH, 'LineStyle', '--', ...
          'Color', 'm');
end

vH = vert_lines([input.delayTimeMs-(0.5*input.rewardWindowWidthMs) input.delayTimeMs input.delayTimeMs+(0.5*input.rewardWindowWidthMs)]);
set(vH, 'Color','k');
grid on
axis tight
ylabel('Number of Trials')
xlabel('Lever Hold Time (ms)')
title('Trial Hold Times Histogram');

%%%%%%%%%%%%%%%%

%% 2 - react time CDF
axH = subplot(subplotSz{:}, 8);
cdfplot([input.reactTimesMs{:}]);
set(gca, 'XLim', [-1000 1000], ...
         'YLim', [0 1]);
hold on;
%This change adds green vertical lines to mark off the reward window, now
%centered around the target time.
vH = vert_lines([-0.5*input.rewardWindowWidthMs 0.5*input.rewardWindowWidthMs]);
set(vH, 'Color','g');
title('Reaction Time CDF');
grid on
xlabel('Time from Reward Window (ms)');

%%%%%%%%%%%%%%%%

%% 3 - React Time PDF
axH = subplot(subplotSz{:}, 6);
nPts = length(input.reactTimesMs);
visIx = reactV<=maxX;
nVisPts = sum(visIx);
if nVisPts > 50
  binWidth = iqr(reactV(visIx))./nVisPts.^(1/3); % robust version of std
  nBins = ceil(2000./binWidth);
else
  nBins = 10;
end
if nBins < 10; nBins=10; end

edges = linspace(-1000, 1000, nBins);
binSize = edges(2)-edges(1);

emptyIx = cellfun(@isempty, input.reactTimesMs);   % see above holdTimesM
if sum(emptyIx) > 0, input.reactTimesMs{emptyIx} = NaN; end

Ns = histc(reactV(successIx), edges);
Nf = histc(reactV(earlyIx), edges);
Nl = histc(reactV(lateIx), edges);
if sum(Ns)+sum(Nf)+sum(Nl) > 0
  bH = bar(edges+binSize/2, [Nf(:),Ns(:), Nl(:)], 'stacked');
  set(bH, 'BarWidth', 1, ...
          'LineStyle', 'none');
  cMap = get(gcf, 'Colormap');
  % flip colors, keep blue on top of red, note flipped in bar.m above
  set(bH(1), 'FaceColor', 'c');
  set(bH(2), 'FaceColor', 'g'); 
  set(bH(3), 'FaceColor', 'm');
end
vH = vert_lines([-0.5*input.rewardWindowWidthMs 0.5*input.rewardWindowWidthMs]);
set(vH, 'Color','k');

hold on;
yLim = get(gca, 'YLim');
plot([0 0], yLim, 'k');
grid on
set(gca, 'XLim', [-1010 1010]);
title('Reaction Times Histogram');
ylabel('Number of Trials');
xlabel('Time from Reward Window (ms)')
%%%%%%%%%%%%%%%%

%% 4 - smoothed perf curve
axH = subplot(subplotSz{:},2:3);
hold on;
%plot(smooth(double(successIx), ceil(nTrial/10), smoothType));
lH = plot(smooth(double(successIx), nTrial, smoothType));
set(lH, 'Color', 'r', ...
        'LineWidth', 3);
lH2 = plot(smooth(double(successIx), 15, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(smooth(double(lateIx), 15, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2, ...
         'LineStyle', '-.');

lH4 = plot(smooth(double(earlyIx), 15, smoothType));
  set(lH4, 'Color', 0.8*[0 1 1], ...
           'LineWidth', 2, ...
           'LineStyle', '-.');
  anystack(lH4, 'bottom');

lH5 = plot(smooth(double(noreleaseIx), 15, smoothType));
  set(lH5, 'Color', orange, ...
           'LineWidth', 2, ...
           'LineStyle', '-.');
  anystack(lH5, 'bottom');

anystack(lH3, 'bottom');


title('Outcomes Plotted by Trials')
ylabel('Percent of Trials');
set(gca, 'YLim', [0 1]);

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
    
      

%%%%%%%%%%%%%%%%

%% 6 - trial length plot
axH = subplot(subplotSz{:}, 7);
hold on;
holdStarts = double(cellvect2mat_padded(input.holdStartsMs));
hSDiffsSec = diff(holdStarts)/1000;

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
  for iT=1:nTrs
    tJ = input.juiceTimesMsCell{iT};
    if isempty(tJ), tJ = 0; end
    firstRewMs(iT) = tJ(1);
    totalRewMs(iT) = sum(tJ);
  end

  rateNSmooth = 45;
  avgRatePerMin = smooth(trPerMin, rateNSmooth, smoothType);
  avgCorrPerMin = smooth(trPerMin.*successIx(2:end), rateNSmooth*1.5, smoothType);
  avgRewPerMin = smooth(trPerMin.*totalRewMs(2:end), rateNSmooth*1.5, smoothType);
  

  [axesH pH1 pH1a]=plotyy(xs,hSCapped, xs, avgRatePerMin);
  set(pH1, 'LineStyle', 'none', ...
           'Marker', '.');
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

ylabel('Trial Start Time (s)');
xLim = trXLim;
set(axesH(~isnan(axesH)), 'XLim', xLim);
lH = plot(xLim, 20*[1 1], '--k');

set(axesH(1),'YLim', [0 121], ...
             'YTick', [0 40 80 120], ...
             'YTickLabelMode', 'auto', ...
             'YColor', 'k');

if ~isnan(axesH(2)) 
  ylabel(axesH(2), 'Trials/min; Avg Reward (ms/sec)');
  yLim2 = get(axesH(2), 'YLim');
  yLim2 = [0 max(yLim2(2), 6)];
  set(axesH(2), 'YLim', yLim2, ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto');
end
  
nDiffs = length(hSDiffsSec);
fN = max(1, nDiffs-5);  % if first 6 trials, start at 1
title(sprintf('Last 6 (sec): %s', mat2str(round(hSDiffsSec(fN:end)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6 - Hold and React Times Over Trials
axH = subplot(subplotSz{:}, 9);
hold on;

% do smooth on only corrects and earlies, then plot by true trial number
desIx = successIx|earlyIx|lateIx;
xVals = find(desIx);
v1 = smooth(holdV(desIx), 25, 'rloess');
v2 = smooth(holdV(desIx), 250, 'rlowess');
desYIx = successIx;
xYVals = find(desYIx);
vy1 = smooth(reactV(successIx), 25, 'rloess');
vy2 = smooth(reactV(successIx), 250, 'rlowess');

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
  ylabel('Hold time (ms) - Excluding NoReleases');
end

% 2nd axes
if ~isempty(vy1) && ~isempty(vy2)
    hyH(1) = plot(axH(2), xYVals, vy1);
    hyH(2) = plot(axH(2), xYVals, vy2);
    
    c2 = 'b';
    set(hyH, 'Color', c2);
    set(hyH(2), 'LineWidth', 2);
    
    set(axH(2),... 'YLim', [0 max(reactV)], ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto', ...
                'YColor', c2)
    ylabel(axH(2), 'React time (ms) - Corrects');
end

title('Mean Reaction (blue) & Hold (black) Times');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic code below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List changed text
axH = subplot(subplotSz{:}, 4);
hold on;
set(axH, 'Visible', 'off');
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
    block2TrPer80Changed = false;
    for iR = 1:nRows
      tVarName = tChangedList{iR,1};
      tChangedTo = tChangedList{iR,2};
      % special case odds changes
      if strfind(tVarName, 'trPer80Level')
        trPer80Changed = true;  % and iterate; summarize below
      elseif strfind(tVarName, 'block2TrPer80Level')
        block2TrPer80Changed = true; % and iterate; summarize below
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
    if block2TrPer80Changed && tDesN > 1
      tStr = mat2str(input.block2TrPer80History{tDesN});
      changedStrOut{end+1} = sprintf('Tr %3d - Block2TrPer80LevelN: -> %s', ...
                                     tDesN, tStr);
    end
  end
  input.changedStr = changedStrOut;

  text(0,1,input.changedStr, ...
       'HorizontalAlignment', 'left', ...
       'VerticalAlignment', 'top');
end


%% Add a save button
if ~exist(consts.behavPdfPath, 'file')
  mkdir(consts.behavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  consts.behavPdfPath, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfunctions

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfig_print(epParamsIn{:});

function outM = subCell2PadVect(cellV, padChar, nMaxTrials)
if nargin < 2, padChar = NaN; end
if nargin < 3, nMaxTrials = length(cellV); end
isEIx = cellfun(@isempty, cellV);
outM = repmat(padChar, [1 nMaxTrials]);
outM(~isEIx) = cell2mat(cellV(~isEIx));

%%%%%%%%%%%%%%%%

function ebH = subManualErrorbar(lH, upperLength, lowerLength)

if nargin < 3, lowerLength = upperLength; end
xv = get(lH, 'XData');
yv = get(lH, 'YData');
xvals = cat(1, xv, xv, xv);
yvals = cat(1, yv-lowerLength(:)', yv+upperLength(:)', yv*NaN);
xvals = reshape(xvals, [1 prod(size(xvals))]);
yvals = reshape(yvals, [1 prod(size(yvals))]);
ebH = plot(xvals, yvals);
set(ebH, 'Color', get(lH, 'Color'));
return

%%%%%%%%%%%%%%%%

function outStr = subPpGrating(input, doPpBlock2Grating)
if nargin < 2, doPpBlock2Grating = false; end

% assemble vars- block2 or not
tV = struct('gratingWidthDeg', []);
if doPpBlock2Grating
  tV.gratingWidthDeg = input.block2GratingWidthDeg;
  tV.gratingHeightDeg = input.block2GratingHeightDeg;
  tV.gratingAzimuthDeg = input.block2GratingAzimuthDeg;
  tV.gratingElevationDeg = input.block2GratingElevationDeg;
  tV.gratingSpatialFreqCPD = input.block2GratingSpatialFreqCPD;
  tV.gratingDurationMs = input.block2GratingDurationMs;
  tV.gratingSpeedDPS = input.block2GratingSpeedDPS;
else
  tV.gratingWidthDeg = input.gratingWidthDeg;
  tV.gratingHeightDeg = input.gratingHeightDeg;
  tV.gratingAzimuthDeg = input.gratingAzimuthDeg;
  tV.gratingElevationDeg = input.gratingElevationDeg;
  tV.gratingSpatialFreqCPD = input.gratingSpatialFreqCPD;
  tV.gratingDurationMs = input.gratingDurationMs;
  tV.gratingSpeedDPS = input.gratingSpeedDPS;
end

outStr = sprintf('%gx%gdeg, at (%g,%g), %gcpd, %dms', ...
                 roundn(double(tV.gratingWidthDeg),-2), ...
                 roundn(double(tV.gratingHeightDeg),-2), ...
                 tV.gratingAzimuthDeg, ...
                 tV.gratingElevationDeg, ...
                 tV.gratingSpatialFreqCPD, ...
                 tV.gratingDurationMs);


if tV.gratingSpeedDPS ~= 0 
  outStr = strcat(outStr, sprintf(' %d deg/s', tV.gratingSpeedDPS));
end

%%%%%%%%%%%%%%%%

function outStr = subPpLaser(input)

if input.laserRampDoExpRamp == 0, 
  rampShapeStr = 'lin'; 
else
  rampShapeStr = 'exp';
end

if input.laserDoLinearRamp
  outStr = [ sprintf('ramp %d (%s)+ const %d ms', ...
                         input.laserRampLengthMs, ...
                         rampShapeStr, ...
                         input.laserRampExtraConstantLengthMs) ];
elseif input.laserDoPulseTrain
  outStr = [ sprintf('train, pulse %d, period %d, dur %d', ...
                         input.laserPulseLengthMs, input.laserPulsePeriodMs, ...
                         input.laserTrainLengthMs) ];
end


%%%%%%%%%%%%%%%%

function outStr = subPpTrialLaser(input, block2Num)

tV(1).trialLaserPowerMw = [];
if block2Num == 1
  prefix = 'trialLaser';
elseif block2Num == 2
  prefix = 'block2TrialLaser';
else
  error('bug 1a');
end
tV.trialLaserPowerMw = input.([prefix 'PowerMw']);
tV.trialLaserOnTimeMs = input.([prefix 'OnTimeMs']);
tV.trialLaserOffTimeMs = input.([prefix 'OffTimeMs']);

outStr = sprintf('%gmW, on %d off %d', ...
                 chop(tV.trialLaserPowerMw, 2), ...
                 chop(tV.trialLaserOnTimeMs, 2), ...
                 chop(tV.trialLaserOffTimeMs, 2));



