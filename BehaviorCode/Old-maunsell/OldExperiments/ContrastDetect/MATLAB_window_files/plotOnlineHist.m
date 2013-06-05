function plotOnlineHist(data_struct, input)

figNum = 4;

wasFigure = ishandle(figNum);

%% draw figure
figH = figure(figNum);
clf;

% set up figure position and style if this is the first launch
if ~wasFigure
  set(figH, 'WindowStyle', 'normal');
  switch hostname
   case 'MaunsellMouse1'
    figPos = [1111 338 806 768];
   otherwise
    figPos = [930 270 780 750];
  end
  set(figH, 'Position', figPos);
end

spSz = {3,3};

% some data processing
nPts = length(input.holdTimesMs);
nTrial = length(input.trialOutcomeCell);
% sometimes on client restart we have empty elements here; pad with
% NaN and hope they match up; I should figure out why they are missing
emptyIx = cellfun(@isempty, input.holdTimesMs); 
if sum(emptyIx) > 0, input.holdTimesMs{emptyIx} = NaN; end
holdV = [input.holdTimesMs{:}];

successIx = strcmp(input.trialOutcomeCell, 'success');
failureIx = strcmp(input.trialOutcomeCell, 'failure');
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');

%% 1 - total hold time histogram
axH = subplot(spSz{:},1);
maxFail = max(holdV(find(failureIx)));
if maxFail > 2000
  maxX = 3000;
elseif maxFail > 2500
  maxX = 2500;
else
  maxX = 2000;
end
visIx = holdV<=maxX;
nVisPts = sum(visIx);
if nVisPts > 50
  binWidth = iqr(holdV(visIx))./nVisPts.^(1/3); % robust version of std
  nBins = ceil(maxX./binWidth);
else
  nBins = 10;
end
edges = linspace(0, maxX, nBins);
Ns = histc(holdV(find(successIx)), edges);
Nf = histc(holdV(find(failureIx)), edges);
if sum(Ns)+sum(Nf) > 0
  bH = bar(edges, [Ns(:),Nf(:)], 'stacked');
  set(bH, 'BarWidth', 1, ...
          'LineStyle', 'none');
end
hold on;
xLim = [0 maxX+50];
set(gca, 'XLim', xLim);
yLim = get(gca, 'YLim');

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
title(sprintf('median hold %4.1fms, react %4.1fms', ...
              median([input.holdTimesMs{:}]), mean([input.reactTimesMs{:}])));
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
if sum(emptyIx) > 0, [input.reactTimesMs{emptyIx}] = deal(NaN); end
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
ylabel('pct correct');
set(gca, 'YLim', [0 1]);


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
xLim = get(gca, 'XLim');
lH = plot(xLim, 20*[1 1], '--k');
set(gca,'YLim', [0 121]);
%%%%%%%%%%%%%%%%

%% 5 - trials per unit time
axH = subplot(spSz{:},5);
hold on;
holdStarts = [input.holdStartsMs{:}];
%edges = linspace(min(holdStarts), max(holdStarts), 100);
%trsPerMin = histc(edges, holdStarts);
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

set(gca,'YLim', [0 61]);

xlabel('trial number');
ylabel('Trials per min');

nDiffs = length(hSDiffsSec);
fN = max(1, nDiffs-5);  % if first 6 trials, start at 1
title(sprintf('Last 6 (sec): %s', mat2str(round(hSDiffsSec(fN:end)))));
%%%%%%%%%%%%%%%%

% 6 - hold times over time
axH = subplot(spSz{:}, 8);
hold on;
hH(1) = plot(smooth(holdV, 50, 'loess'));
hH(2) = plot(smooth(holdV, 250, 'loess'));
set(hH(2), 'Color', 'k', ...
           'LineWidth', 3);
hRange = [0 prctile(holdV, 95)];
set(gca, 'YLim', hRange);
xlabel('trial number');
ylabel('Hold time (ms)');
if nDiffs>0
  totalElapsedS = (input.holdStartsMs{end} - ...
                   input.holdStartsMs{1})/1000;

  %totalRewMs = sum([input.totalRewardTimesMs{successIx}]);
  totalRewMs = sum(cat(2,input.juiceTimesMsCell{:}));
  title(sprintf('Total elapsed: %dmin; reward %.1fsec', ...
                round(totalElapsedS/60), totalRewMs/1000));
end
%%%%%%%%%%%%%%%%

% 7 - contrast performance
axH = subplot(spSz{:}, 3);

% get contrast data
ds = load('~/Library/Application Support/MWorks/ContrastDetectGeneratedImages/info.mat');
contrastP = struct('correct', repmat(NaN, [1 ds.nContrasts]), ...
                   'total', repmat(NaN, [1 ds.nContrasts]));
for iC = 1:ds.nContrasts
  imgIx = cellfun(@(x)isequal(x,iC-1), input.imageNumber);  % 0 origin in XML
  contrastP.total(iC) = sum(imgIx);
  contrastP.correct(iC) = sum(imgIx & successIx);
end
contrastP.fractCorrect = contrastP.correct./contrastP.total;
contrastP.sem = (contrastP.fractCorrect.*contrastP.fractCorrect) ./ sqrt(contrastP.total);
% save contrast data
input.contrastP = contrastP; 
input.contrastLevelsPct = ds.contrastLevelsPct;

% draw figure
plot(input.contrastLevelsPct, contrastP.fractCorrect*100, 'x-');
% error bars
hold on;
b1H = plot(ds.contrastLevelsPct, (contrastP.fractCorrect+contrastP.sem)*100);
b2H = plot(ds.contrastLevelsPct, (contrastP.fractCorrect-contrastP.sem)*100);
set([b1H;b2H], 'Color', 0.8*[1 1 1]);
xlabel('Contrast (%)');
ylabel('Correct (%)');

% $$$ set(gca, 'XLim', [0.1 100], ...
% $$$          'XTick', [0.1 1 10 100], ...
% $$$          'XTickLabel', {'0.1', '1', '10', '100'});
%xtickloc = get(gca, 'XTick')
%
set(gca, 'YLim', [0 100]);
set(gca, 'XScale', 'log');

% autocompute xlimits and ticks
xLim = [0.5 1.5] .* [min(ds.contrastLevelsPct), max(ds.contrastLevelsPct)]
set(gca, 'XLim', xLim);
xTAll = [ (1:9)/100, (1:9)/10, (1:9), (1:9)*10, 100];
desXTIx = xTAll >= xLim(1) & xTAll <= xLim(2);
xticks = xTAll(desXTIx);

grid on;
set(gca, 'XTick', xticks)

% labels
xtLabIx = xticks == 0.01 | xticks == 0.1 | xticks == 1 | xticks == 10 ...
         | xticks == 100;
xtLab = xticks(xtLabIx);
labStrs = cell(length(xticks), 1);
labels = cellstr(num2str(xtLab(:)));
labStrs(xtLabIx) = labels;
set(gca, 'XTickLabel', labStrs);

%%%%%%%%%%%%%%%%

%% supertitle
nCorr = sum(successIx);
nFail = sum(failureIx);
nIg = sum(ignoreIx);
suptitle(sprintf('%03d: Correct %d (%4.1f%%), early %d (%4.1f%%), ign %d (%4.1f%%) total %d', ...
                 input.subjectNum, ...
                 nCorr, nCorr/nTrial*100, nFail, nFail/nTrial*100, ...
                 nIg, nIg/nTrial*100, nTrial));

%% Add a save button

outName = sprintf('/Users/histed/behav-output/%s-behav-i%03d.pdf', ...
                  datestr(now, 'yymmdd'), input.subjectNum);
epParams = { figNum, outName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
bH = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF figure : %s', outName), ...
               'Units', 'pixels', ...
               'Position', [5 5 450 20], ...
               'Callback', { @saveButtonCb, epParams });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfunctions

function saveButtonCb(hObject, eventdata, epParamsIn) 
exportfigPrint(epParamsIn{:});

%%%%%%%%%%%%%%%%

