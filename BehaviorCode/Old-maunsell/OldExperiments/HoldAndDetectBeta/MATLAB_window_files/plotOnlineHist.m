function input = plotOnlineHist(data_struct, input)

figNum = 4;
name = 'HoldAndDetectBeta';
spSz = {4,3};

%% draw figure
figH = figure(figNum);
set(figH, 'ToolBar', 'none');
clf;

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
% sometimes on client restart we have empty elements here; pad with
% NaN and hope they match up; I should figure out why they are missing
emptyIx = cellfun(@isempty, input.holdTimesMs); 
if sum(emptyIx) > 0, input.holdTimesMs{emptyIx} = NaN; end
holdV = [input.holdTimesMs{:}];
t1 = input.reqHoldTimeMs;
[t1{cellfun(@isempty, t1)}] = deal(NaN);
reqHoldV = [t1{:}];

successIx = strcmp(input.trialOutcomeCell, 'success');
failureIx = strcmp(input.trialOutcomeCell, 'failure');
earlyIx = failureIx;
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
nCorr = sum(successIx);
nFail = sum(failureIx);
nIg = sum(ignoreIx);
holdStarts = [input.holdStartsMs{:}];
nTrials = length(input.trialOutcomeCell);


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
	axH = subplot(spSz{:}, 1);						% default axes are 0 to 1
	set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])
        numTrials = nTrial;
	text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 14);
	text(0.00, 1.0, {'Subject:', 'Time working:', 'Reward volume:'});
	text(0.60, 1.0, {sprintf('%d', input.subjectNum), ...
		sprintf('%.0d m', ...
		round((input.holdStartsMs{end} - input.holdStartsMs{1})/60000)),...
		sprintf('%.1f s', ...
		sum(cat(2,input.juiceTimesMsCell{:})) / 1000)});
	text(0.00, 0.6, {'Trials:', 'Correct:', 'Early:', 'Failed:'});
	text(0.40, 0.6, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
				sprintf('%d', nFail), sprintf('%d', nIg)});
	text(0.60, 0.6, {' ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
				sprintf('%.0f%%', nFail / numTrials * 100.0), ...
				sprintf('%.0f%%', nIg / numTrials * ...
                                        100.0)});

        text(0.0, 0.32, sprintf(['Hold, react median:']));

        text(0.65, 0.32, sprintf(['%5.1f, %5.1f ms'], ...
                               median([input.holdTimesMs{:}]), ...
                               median([input.reactTimesMs{:}])), ...
             'HorizontalAlignment', 'left');
        
        text(0.0, 0.2, sprintf(['Hold (f+r, tf): \t%3d+%3d,%3d ms\n', ...
                                'Timeouts (e,m):\t%4.1f, %4.1f s\n', ...
                                'React:\t%5.2f s\n'], ...
                                input.fixedReqHoldTimeMs, ...
                                input.randReqHoldMaxMs, ...
                                input.tooFastTimeMs, ...
                                input.earlyTimeoutMs/1000, ...
                                input.missedTimeoutMs/1000, ...
                                input.reactTimeMs/1000), ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Hold time histogram

axH = subplot(spSz{:}, 3);
maxFail = max(holdV(find(failureIx|successIx)));
if maxFail >= 2500
  maxX = max(holdV(successIx))+200;
elseif maxFail >= 2000
  maxX = 2500;
else
  maxX = 2000;
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
maxXHoldAxis = maxX+300;
xLim = [0 maxXHoldAxis];
set(gca, 'XLim', xLim);
yLim = get(gca, 'YLim');

if ~isempty(Nf)
  pH = plot(edges, Nf);
  set(pH, 'LineStyle', '--', ...
          'Color', 'r');
end

title('total hold times');


if ~isempty(input.tooFastTimeMs)
  yLim = get(gca, 'YLim');
  plot(input.tooFastTimeMs * [1 1], yLim, 'k--');
end
%%%%%%%%%%%%%%%%

%% 2 - react time CDF
axH = subplot(spSz{:},4);
cdfplot([input.reactTimesMs{:}]);
set(gca, 'XLim', [-1000 1000], ...
         'YLim', [0 1]);
hold on;
vH = vert_lines(0:200:1000);
set(vH, 'LineStyle', ':', 'Color', 0.5*[1 1 1]);
title('');
xlabel('');

%%%%%%%%%%%%%%%%

%% 3 - react time PDF
axH = subplot(spSz{:},7);
nPts = length(input.reactTimesMs);
reactV = [input.reactTimesMs{:}];
visIx = reactV<=maxX;
nVisPts = sum(visIx);
if nVisPts > 50
  binWidth = iqr(reactV(visIx))./nVisPts.^(1/3); % robust version of std
  nBins = ceil(2000./binWidth);
else
  nBins = 10;
end
edges = linspace(-1000, 1000, nBins);
binSize = edges(2)-edges(1);

emptyIx = cellfun(@isempty, input.reactTimesMs);   % see above holdTimesM
if sum(emptyIx) > 0, input.reactTimesMs{emptyIx} = NaN; end
rV = [input.reactTimesMs{:}];

Ns = histc(rV(successIx), edges);
Nf = histc(rV(failureIx), edges);
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
title('reaction times');
%%%%%%%%%%%%%%%%

%% 4 - smoothed perf curve
axH = subplot(spSz{:},2);
hold on;
plot(smooth(double(successIx), ceil(nTrial/10), 'lowess'));
lH = plot(smooth(double(successIx), nTrial, 'lowess'));
set(lH, 'Color', 'r', ...
        'LineWidth', 3);
lH2 = plot(smooth(double(successIx), 100, 'lowess'));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);

lH3 = plot(smooth(double(ignoreIx), 100, 'lowess'));
set(lH3, 'Color', 'm', ...
         'LineWidth', 2, ...
         'LineStyle', ':');
lH4 = plot(smooth(double(failureIx), 100, 'lowess'));
set(lH4, 'Color', 'c', ...
         'LineWidth', 2, ...
         'LineStyle', ':');
anystack(lH4, 'bottom');
anystack(lH3, 'bottom');


ylabel('pct correct');
set(gca, 'YLim', [0 1]);

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
    
      

%%%%%%%%%%%%%%%%

%% 6 - trial length plot
axH = subplot(spSz{:},6);
hold on;
holdStarts = [input.holdStartsMs{:}];
%pH=semilogy(diff(holdStarts)/1000);
hSDiffsSec = diff(holdStarts)/1000;
% make outliers a fixed value
largeIx = hSDiffsSec >= 120;
hSDiffsSec(largeIx) = 120;
xs = 1:length(hSDiffsSec);
pH1=plot(xs,hSDiffsSec);
if sum(largeIx) > 0
  pH2 = plot(xs(largeIx),hSDiffsSec(largeIx),'r.');  % outliers
end
set(pH1, 'LineStyle', 'none', ...
        'Marker', 'x');
pH2 = plot(smooth(hSDiffsSec, 5, 'lowess'), 'r');

%plot(diff(holdStarts)/1000, 'x');
ylabel('trial start time diff (s)');
xLim = trXLim;
set(gca, 'XLim', xLim);
lH = plot(xLim, 20*[1 1], '--k');
set(gca,'YLim', [0 121]);

nDiffs = length(hSDiffsSec);
fN = max(1, nDiffs-5);  % if first 6 trials, start at 1
title(sprintf('Last 6 (sec): %s', mat2str(round(hSDiffsSec(fN:end)))));

%%%%%%%%%%%%%%%%

%% 5 cumulative time elapsed
axH = subplot(spSz{:},5);
hold on;

xs = 1:length(hSDiffsSec);
pH1 = plot(xs, cumsum(hSDiffsSec)./60, '.-');
maxMin = sum(hSDiffsSec)./60;

ylabel('Time working (min)');
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
axH = subplot(spSz{:}, 8);
hold on;

% do smooth on only corrects and earlies, then plot by true trial number
desIx = successIx|earlyIx;
xVals = find(desIx);
v1 = smooth(holdV(desIx), 25, 'loess');
v2 = smooth(holdV(desIx), 250, 'lowess');

if ~isempty(v1) && ~isempty(v2)
  hH(1) = plot(xVals, v1);
  hH(2) = plot(xVals, v2);
  set(hH(2), 'Color', 'k', ...
           'LineWidth', 3);
end
%hRange = [0 prctile(holdV(desIx), 95)];
set(gca, 'XLim', trXLim);
yLim = get(gca, 'YLim');
if yLim(2) > 0
  yLim(1) = 0;
  set(gca, 'YLim', yLim);
end
ylabel('Hold time (ms)');
if nDiffs>0
  totalElapsedS = (input.holdStartsMs{end} - ...
                   input.holdStartsMs{1})/1000;

  %totalRewMs = sum([input.totalRewardTimesMs{successIx}]);
  totalRewMs = sum(cat(2,input.juiceTimesMsCell{:}));

  title('mean hold over corrects+early');
  
end

% 9 - reward sizes over correct trials
nTrs = length(input.juiceTimesMsCell);
firstSize = repmat(NaN, [1 nTrs]);
totalSize = repmat(NaN, [1 nTrs]);
for iT=1:nTrs
  tJ = input.juiceTimesMsCell{iT};
  if isempty(tJ), tJ = 0; end
  firstSize(iT) = tJ(1);
  totalSize(iT) = sum(tJ);
end

% 10 -  number of responses during react window
axH = subplot(spSz{:}, 10);
winLen = 175;
winStart = 250;
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
title('testing')

%xLim = [0 nTrials];
xLim = trXLim;
set(gca, 'YLim', [0 1], ...
         'XLim', xLim);
plot(xLim, 0.5*[1 1], '--');

axH = subplot(spSz{:}, 4);
hold on;
v2H = vert_lines([winStart winStart+winLen]);
set(v2H, 'Color', 'g', 'LineStyle', '--');

%%%%%%%%%%%%%%%%

%% 11 - mean reward over correct and early trials
 
axH = subplot(spSz{:}, 11);
hold on;

desIx = successIx | earlyIx;
xv = find(desIx);
sH1 = plot(xv, smooth1(firstSize(desIx), 'gauss', 1, 11), 'b');
sH2 = plot(xv, smooth1(totalSize(desIx), 'gauss', 1, 11), 'r');
sH3 = plot(xv, smooth1(firstSize(desIx), 'gauss', 2, 51), 'b');
sH4 = plot(xv, smooth1(totalSize(desIx), 'gauss', 2, 51), 'r');
set([sH3;sH4], 'LineWidth', 3);
ylabel('Reward size (ms)');
title('reward over corrects+early');
set(gca, 'XLim', trXLim);
xlabel('Trial number');

%%%%%%%%%%%%%%%%

%% 12: error rate with hold time?
axH = subplot(spSz{:}, 12);
maxReqX = max(input.fixedReqHoldTimeMs)+ ...
          max(input.randReqHoldMaxMs);
if nVisPts > 50
  binWidth = iqr(reqHoldV(visIx)) ./ nVisPts.^(1/3);	% robust version of std
  nBins = ceil(maxReqX ./ binWidth);
else
  nBins = 10;
end
edges = linspace(0,maxReqX,nBins);

zVec = edges*0;
if sum(successIx) == 0
  N1 = zVec;
else
  N1 = histc(reqHoldV(successIx), edges); 
end
if sum(earlyIx) == 0
  N2 = zVec;
else
  N2 = histc(reqHoldV(earlyIx), edges);
end
if sum(ignoreIx) == 0
  N3 = zVec;
else
  N3 = histc(reqHoldV(ignoreIx), edges);
end

totN = N1+N2+N3;
totN(totN==0) = Inf;

hold on;
clear('pH')
pH(1) = plot(edges, N1./totN*100, 'k');
pH(2) = plot(edges, N2./totN*100, 'c--');
pH(3) = plot(edges, N3./totN*100, 'm--');
set(pH, 'LineWidth', 2)

set(gca, 'XLim', [0 maxReqX*1.1]);
set(gca, 'YLim', [0 100]);
ylabel('% of trials');
title('trial outcome by hold time');

%%%%%%%%%%%%%%%%
%% testing plots

doTestingFig = false;

if doTestingFig
  figTestH = figure(5);
  clf;

  %% trials per unit time
  testSpSz={2,2};
  axH = subplot(testSpSz{:},1);
  hold on;
  holdStarts = [input.holdStartsMs{:}];
  hSDiffsSec = diff(holdStarts)/1000;
  trsPerSec = 1./hSDiffsSec;
  trsPerMin = trsPerSec*60;
  xs = 1:length(hSDiffsSec);

  h1 = plot(xs, trsPerMin);
  set(h1, 'LineStyle', 'none', ...
          'Marker', 'x', ....
          'Color', 'b');

  hold on;
  pH2 = plot(smooth(trsPerMin, 50, 'loess'));
  set(pH2, 'Color', 'k', ...
           'LineWidth', 3);

  % line
  xLim = get(gca, 'XLim');
  lH = plot(xLim, 20*[1 1], '--k');

  ylabel('Trials per min');

  %% elapsed
  axH = subplot(testSpSz{:},2);
  hold on;

  holdStarts = [input.holdStartsMs{:}];
  hSDiffsSec = diff(holdStarts)/1000;
  xs = 1:length(hSDiffsSec);
  pH1=plot(xs,hSDiffsSec);
  set(pH1, 'LineStyle', 'none', ...
           'Marker', 'x');
  pH2 = plot(smooth(hSDiffsSec, 5, 'lowess'), 'r');

  ylabel('trial start time diff (s)');
  xLim = trXLim;
  set(gca, 'XLim', xLim);
  set(gca,'YLim', [0 20]);

  %% end testing plot
  suptitle2('testing');
  figH = figure(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic code below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List changed text
axH = subplot(spSz{:}, 9);
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
    for iR = 1:nRows
      tVarName = tChangedList{iR,1};
      tChangedTo = tChangedList{iR,2};

      changedStrOut{end+1} = sprintf('Tr %3d - %s: -> %g', ...
                                     tDesN, tVarName, tChangedTo);
    end
  end
  input.changedStr = changedStrOut;

  text(0,1,input.changedStr, ...
       'HorizontalAlignment', 'left', ...
       'VerticalAlignment', 'top');
end

%%%%%%%%%%%%%%%%


%% Add a save button

outDir = '/Users/holdanddetect/Documents/MWorks/BehavOutputPdfs';
if ~exist(outDir, 'file')
  mkdir(outDir)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  outDir, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfunctions

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfigPrint(epParamsIn{:});

%%%%%%%%%%%%%%%%

