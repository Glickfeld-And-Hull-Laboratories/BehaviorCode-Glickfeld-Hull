function input = plotRewardDelay(data_struct, input)

% essential consts
figNum = 4;
name = 'RewardDelay';
consts = reward_delay_constants;
subplotSz = {4,3};

smoothType = 'lowess';

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
nPts = length(input.holdTimesMs);
nTrial = length(input.trialOutcomeCell);
reqHoldV = celleqel2mat_padded(input.tTotalReqHoldTimeMs);
successIx = strcmp(input.trialOutcomeCell, 'success');
failureIx = strcmp(input.trialOutcomeCell, 'failure');
earlyIx = failureIx;
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
missIx = ignoreIx;
nCorr = sum(successIx);
nFail = sum(failureIx);
nIg = sum(ignoreIx);
holdStarts = [input.holdStartsMs{:}];
nTrials = length(input.trialOutcomeCell);
tGotLongBonusThisTrialV = celleqel2mat_padded(input.tGotLongBonusThisTrial);

reactV = celleqel2mat_padded(input.reactTimesMs);
holdV = celleqel2mat_padded(input.holdTimesMs);  % sometimes on client restart have empty elements here; should figure out why
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;
tTrialN = input.trialSinceReset;

%fprintf(1,'%d %d\n', sum(block2Tr1Ix&~earlyIx), sum(block2Tr2Ix&~earlyIx));  %% debug # of block2 trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Values

if length(holdStarts) > 2
	holdTimeMS = input.holdTimesMs{end};
	switch input.trialOutcomeCell{end}
	case 'success'
		outcomeString = 'correct';
	case 'failure'
		outcomeString = 'early';
	case 'ignore'
		outcomeString = 'failed';
	end

        axH = subplot(subplotSz{:}, 1);						% default axes are 0 to 1
        
	set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

        numTrials = nTrial;
        
        set(gcf, 'Visible', 'off'); % hide figure during text
                                    % drawing - kludge

        text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 16);

        elMin = round((now - datenum(input.startDateVec)) * 24*60);
        startStr = datestr(input.startDateVec, 'HH:MM');
        text(0.00, 1.05, {'Subject:', 'Start time + elapsed:', 'Reward vol (m\pm std):'});
	text(0.60, 1.05, ...
             { sprintf('%2d', input.subjectNum), ...
               sprintf('%s + %2dm', ...
                       startStr, elMin), ...
               sprintf('%.1f s     \t(%g ms\\pm %g ms)', ...
                       nansum(juiceTimesMsV./1000), ...
                       chop(nanmean(juiceTimesMsV),2), ...
                       chop(nanstd(juiceTimesMsV),2)), ...
             });
	t2H(1) = text(0.00, 0.9, {'Trials:', 'Correct:', 'Early:', 'Failed:'});
	t2H(2) = text(0.35, 0.9, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
				sprintf('%d', nFail), sprintf('%d', nIg)});
	t2H(3) = text(0.54, 0.9, {' ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
				sprintf('%.0f%%', nFail / numTrials * 100.0), ...
				sprintf('%.0f%%', nIg / numTrials * ...
                                        100.0)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');

        
        if input.doLongBonus
          lbStr = sprintf('Long bonus:\n    t +%d, rew +%d ms', ...
                          input.longBonusExtraHoldTimeMs, ...
                          round(input.longBonusExtraRewardUs/ ...
                                1000));
        else
          lbStr = 'Long bonus: off';
        end

        tStr = sprintf(['Hold, react median:\n', ...
                        '   %5.1f ms, %5.1f ms\n', ...
                        '\n', ...
                        '%s'], ...
                       median(holdV), ...
                       median(reactV(successIx)), ...
                       lbStr);
               
        text(0.8, 0.9, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left');
        
        if input.doGeomHoldDist
          randTypeStr = 'geo';
          geoExtraStr = sprintf('(geo m%d)', input.geomHoldMeanMs);
        else 
          randTypeStr = 'Uniform';
          geoExtraStr = '';
        end

        
        tStr = sprintf( ['Hold (f+r/%s, tf): \t%3d+%3d,%3d ms %s\n', ...
                         'Timeouts (e,m):\t%4.1f, %4.1f s\n', ...
                         'React:\t%5.2f s;   ITI %d +%d ms\n', ...
                         'Reward min, max: %3.0f, %3.0f ms   %s\n'], ...
                        randTypeStr, ...
                        input.fixedReqHoldTimeMs, ...
                        input.randReqHoldMaxMs, ...
                        input.tooFastTimeMs, ...
                        geoExtraStr, ...
                        input.earlyTimeoutMs/1000, ...
                        input.missedTimeoutMs/1000, ...
                        input.reactTimeMs/1000, ...
                        input.itiTimeMs, ...
                        input.itiExtraRandTimeMs, ...
                        input.minRewardUs/1000, ...
                        input.maxRewardUs/1000);

        text(0.0, 0.55, tStr, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

        set(gcf, 'Visible', 'on');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Hold time histogram

axH = subplot(subplotSz{:}, 3);
if input.fixedReqHoldTimeMs+input.randReqHoldMaxMs < 1000
  maxX = 2000;
else
  maxFail = double(input.reactTimeMs+input.fixedReqHoldTimeMs+input.randReqHoldMaxMs);
  maxX = ceil((maxFail+45)./500)*500;  % round up to nearest 500 ms.
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
Nf = histc(holdV(find(failureIx)), edges);
if sum(Ns) + sum(Nf) > 0
  bH = bar(edges, [Ns(:),Nf(:)], 'stacked');
  set(bH, 'BarWidth', 1, 'LineStyle', 'none');
end
hold on;
xLim = [0 maxX];
set(gca, 'XLim', xLim);
yLim = get(gca, 'YLim');
if (length(get(gca, 'XTick')) > 4)
  xT = (0:(round(maxX/1000)))*1000;
  set(gca, 'XTick', xT);
end


if ~isempty(Nf)
  pH = plot(edges, Nf);
  set(pH, 'LineStyle', '--', ...
          'Color', 'r');
end

title('Total Hold Times');
if ~isempty(input.tooFastTimeMs)
  yLim = get(gca, 'YLim');
  plot(double(input.tooFastTimeMs) * [1 1], yLim, 'k--');
end
%%%%%%%%%%%%%%%%

%% 2 - react time CDF
axH = subplot(subplotSz{:},4);
cdfplot([input.reactTimesMs{:}]);
set(gca, 'XLim', [-1000 1000], ...
         'YLim', [0 1]);
hold on;
vH = vert_lines(0:200:1000);
set(vH, 'LineStyle', ':', 'Color', 0.5*[1 1 1]);
title('Reaction Time CDF');
xlabel('Time from Reward Window (ms)');
ylabel('Percent of Trials');

%%%%%%%%%%%%%%%%

%% 3 - React Time PDF
axH = subplot(subplotSz{:},7);
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
Nf = histc(reactV(failureIx), edges);
if sum(Ns)+sum(Nf) > 0
  bH = bar(edges+binSize/2, [Nf(:),Ns(:)], 'stacked');
  set(bH, 'BarWidth', 1, ...
          'LineStyle', 'none');
  cMap = get(gcf, 'Colormap');
  % flip colors, keep blue on top of red, note flipped in bar.m above
  set(bH(1), 'FaceColor', [0.6 0 0]);
  set(bH(2), 'FaceColor', [0 0 0.6]);      
end

hold on;
yLim = get(gca, 'YLim');
plot([0 0], yLim, 'k');

set(gca, 'XLim', [-1010 1010]);
title('Reaction Time PDF');
ylabel('Percent of Trials');
xlabel('Time from Reward Window (ms)')
%%%%%%%%%%%%%%%%

%% 4 - smoothed perf curve
axH = subplot(subplotSz{:},2);
hold on;
plot(smooth(double(successIx), ceil(nTrial/10), smoothType));
lH = plot(smooth(double(successIx), nTrial, smoothType));
set(lH, 'Color', 'r', ...
        'LineWidth', 3);
lH2 = plot(smooth(double(successIx), 100, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(smooth(double(ignoreIx), 100, smoothType));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2, ...
         'LineStyle', '-.');

lH4 = plot(smooth(double(failureIx), 100, smoothType));
  set(lH4, 'Color', 0.8*[0 1 1], ...
           'LineWidth', 2, ...
           'LineStyle', '-.');
  anystack(lH4, 'bottom');

anystack(lH3, 'bottom');

if input.doLongBonus
  %tV = smooth(double(tGotLongBonusThisTrialV(successIx)), 50, smoothType);
  %lH5 = plot(find(successIx), tV);
  tV = smooth(double(tGotLongBonusThisTrialV), 100, smoothType);
  lH5 = plot(tV);
  set(lH5, 'Color', 0.5*[1 1 1], ...
           'LineWidth', 2);
end
title('Outcomes Plotted by Trials')
ylabel('Percent of Trials');
set(gca, 'YLim', [0 1]);

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
    
      

%%%%%%%%%%%%%%%%

%% 6 - trial length plot
axH = subplot(subplotSz{:},6);
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

ylabel('trial start time diff (s)');
xLim = trXLim;
set(axesH(~isnan(axesH)), 'XLim', xLim);
lH = plot(xLim, 20*[1 1], '--k');

set(axesH(1),'YLim', [0 121], ...
             'YTick', [0 40 80 120], ...
             'YTickLabelMode', 'auto', ...
             'YColor', 'k');

if ~isnan(axesH(2)) 
  ylabel(axesH(2), 'Trials/min; avg rew (ms/sec)');
  yLim2 = get(axesH(2), 'YLim');
  yLim2 = [0 max(yLim2(2), 6)];
  set(axesH(2), 'YLim', yLim2, ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto');
end
  
nDiffs = length(hSDiffsSec);
fN = max(1, nDiffs-5);  % if first 6 trials, start at 1
title(sprintf('Last 6 (sec): %s', mat2str(round(hSDiffsSec(fN:end)))));

%%%%%%%%%%%%%%%%

%% 5 cumulative time elapsed
axH = subplot(subplotSz{:},5);
hold on;

hSDiffsRealSec = diff(holdStarts)/1000;
xs = 1:length(hSDiffsRealSec);
pH1 = plot(xs, cumsum(hSDiffsRealSec)./60, '.-');
maxMin = sum(hSDiffsRealSec)./60;

ylabel('Total Training Time (min)');
%set(gca, 'DataAspectRatio', [100 10 1]);
xLim = trXLim;
set(gca, 'XLim', xLim);
% compute data aspect ratio manually: matlab will keep y dim of
% plot box fixed and change x and we want the reverse
yLim = [0 max(xLim(2)/10, maxMin*1.1)];   
set(gca, 'YLim', yLim);
pH2 = plot(xLim, xLim/10, 'k--');

%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%

% 6 - hold times over time
axH = subplot(subplotSz{:}, 8);
hold on;

% do smooth on only corrects and earlies, then plot by true trial number
desIx = successIx|earlyIx;
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
  ylabel('Hold time (ms) - corr+early');
end

% 2nd axes
if ~isempty(vy1) && ~isempty(vy2)
    hyH(1) = plot(axH(2), xYVals, vy1);
    hyH(2) = plot(axH(2), xYVals, vy2);
    
    c2 = 'b';
    set(hyH, 'Color', c2);
    set(hyH(2), 'LineWidth', 2);
    
    set(axH(2), 'YLim', [0 max(reactV)], ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto', ...
                'YColor', c2)
    ylabel(axH(2), 'React time (ms) - corr');
end

title('mean react/hold over time (robust fit)');

% used to have % 9 - reward sizes over correct trials


%% 10 -  Percentage of Responses in the Reward Window (Target)
nStims=1;

if nStims == 1
  axH = subplot(subplotSz{:}, 10);
  winLen = 200;
  winStart = 200;
  winReal = reactV > winStart & reactV <= (winStart+winLen);
  baselineRange = [-winLen*2 winLen];
  nBaseWins = range(baselineRange)./winLen;
  winBaseline = reactV > baselineRange(1) & reactV <= baselineRange(2);
  winBefore = winBaseline ./ nBaseWins;
  
  
  hold on;
  p1H = plot(smooth(double(winReal), ceil(nTrial/10), 'lowess'));
  nSm = 5;
  v1 = smooth(double(winReal), ceil(nTrial./nSm), 'lowess');
  v2 = smooth(double(winBefore), ceil(nTrial./nSm), 'lowess');
  v0 = v1 ./ (v1+v2);
  p2H = plot(v0, 'k');
  set(p2H, 'LineWidth', 3);
  %p2H = plot(smooth(double(winBefore), ceil(nTrial/10), 'lowess'));
  %p3H = plot(smooth(double(winAfter), ceil(nTrial/10), 'lowess'));
  title('trs inside react win')
  
  %xLim = [0 nTrials];
  xLim = trXLim;
  set(gca, 'YLim', [0 1], ...
           'XLim', xLim);
  plot(xLim, 0.5*[1 1], '--');
  
  axH = subplot(subplotSz{:}, 4);
  hold on;
  v2H = vert_lines([winStart winStart+winLen]);
  set(v2H, 'Color', 'g', 'LineStyle', '--');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic code below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List changed text
axH = subplot(subplotSz{:}, 9);
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



