function input = plotDualPressHist(data_struct, input)

% essential consts
figNum = 4;
name = 'DualPress';
cs = exptConstants;
spSz = {4,3};

smoothType = 'lowess';

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

%% color consts
cMap = cmap_unique13;
colBrown = cMap(7,:);
colOrange = [1 0.5 0.2];
colPurple = [0.5 0 0.5];
colYellow = [245 208 76]./256;

% some data processing
nPts = length(input.holdTimesMs);
nTrial = length(input.trialOutcomeCell);

successIx = strcmp(input.trialOutcomeCell, 'success');
earlyIx = strcmp(input.trialOutcomeCell, 'early');
incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
missIx = ignoreIx;
dualReleaseIx= strcmp(input.trialOutcomeCell, 'dualrelease');

nDR = sum(dualReleaseIx);
nCorr = sum(successIx);
nEarly = sum(earlyIx);
nIg = sum(ignoreIx);
nInc = sum(incorrectIx);

nTrials = length(input.trialOutcomeCell);

% make vectors from cell vectors
holdStarts = double(cellvect2mat_padded(input.holdStartsMs));
reactV = double(cellvect2mat_padded(input.reactTimesMs));
holdV = double(cellvect2mat_padded(input.holdTimesMs));  % sometimes on client restart have empty elements here; should figure out why

juiceTimesSumMsV = cellfun(@sum, input.juiceTimesMsCell);  % note MH 121218: if using multiple rewards this code needs to be reworked
juiceTimesSumMsV(juiceTimesSumMsV==0) = NaN;

tTotalRewardTimeMs = celleqel2mat_padded(input.tTotalRewardTimeUs, NaN, 'double')/1000;
tLeftTrial = celleqel2mat_padded(input.tLeftTrial, NaN, 'double');
tFirstReactReleaseIsLeft = celleqel2mat_padded(input.tFirstReactReleaseIsLeft, NaN, 'double');
leverLatencyRPosMs = double(cellvect2mat_padded(input.leverLatencyRPosMs));
tRewardRunningLeftBias = celleqel2mat_padded(input.tRewardRunningLeftBias, NaN, 'double');

tSwitchV = [0 diff(tLeftTrial)~=0];

%% needs to be cleaned up: instead of first limiting to left/right
% e.g. leftOutcomes=toc(leftIx), leftCorr = leftOutcomes =='success', 
% use: leftCorr= successIx & leftIx

%both sides outcomes
leftIx= tLeftTrial == 1;
leftOutcomes=input.trialOutcomeCell(leftIx);
leftTrNs = find(leftIx);
rightOutcomes=input.trialOutcomeCell(~leftIx);
rightTrNs = find(~leftIx);

%left outcomes
nLeftTrs = sum(leftIx);
leftCorr=strcmp(leftOutcomes, 'success');
pLeftCorr = sum(leftCorr)/nLeftTrs;
leftEarly = strcmp(leftOutcomes, 'early');
pLeftEarly = sum(leftEarly)/nLeftTrs;
leftIncorr = strcmp(leftOutcomes, 'incorrect');
pLeftIncorr = sum(leftIncorr)/nLeftTrs;
leftIg = strcmp(leftOutcomes, 'ignore');
pLeftIg = sum(leftIg)/nLeftTrs;
leftDR = strcmp(leftOutcomes, 'dualrelease');
pLeftDR = sum(leftDR)/nLeftTrs;

%right outcomes
nRightTrs = sum(~leftIx);
rightCorr = strcmp(rightOutcomes, 'success');
pRightCorr = sum(rightCorr)/nRightTrs;
rightEarly= strcmp(rightOutcomes, 'early');
pRightEarly = sum(rightEarly)/nRightTrs;
rightIncorr = strcmp(rightOutcomes, 'incorrect');
pRightIncorr = sum(rightIncorr)/nRightTrs;
rightIg = strcmp(rightOutcomes, 'ignore');
pRightIg = sum(rightIg)/nRightTrs;
rightDR = strcmp(rightOutcomes, 'dualrelease');
pRightDR = sum(rightDR)/nRightTrs;

leftReactTimes= reactV(leftIx);
rightReactTimes= reactV(~leftIx);

respMadeIx = (successIx | dualReleaseIx | incorrectIx);
leftECIx= leftIx & (earlyIx | respMadeIx);
leftECTrs= find(leftECIx);
leftHolds= holdV(leftECIx);
leftCIx= leftIx & respMadeIx;
leftCTrs= find(leftCIx);
leftReacts= reactV(leftCIx);


rightECIx= ~leftIx & (earlyIx | respMadeIx);
rightECTrs= find(rightECIx);
rightHolds= holdV(rightECIx);
rightCIx= ~leftIx & respMadeIx;
rightCTrs=find(rightCIx);
rightReacts= reactV(rightCIx);


% Crude Latency Processing for L/R Latency Mean Strings
strLTrIx = (successIx | incorrectIx | dualReleaseIx) & tLeftTrial==1;
strRTrIx = (successIx | incorrectIx | dualReleaseIx) & tLeftTrial==0;
levLatencyV = leverLatencyRPosMs;
leftLatencyMeanMs = nanmean(levLatencyV(strLTrIx));
rightLatencyMeanMs = nanmean(levLatencyV(strRTrIx));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(holdStarts) > 2
    holdTimeMS = input.holdTimesMs{end};
    switch input.trialOutcomeCell{end}
      case 'success'
        outcomeString = 'correct';
      case 'early'
        outcomeString = 'early';
      case 'ignore'
        outcomeString = 'failed';
      case 'incorrect'
        outcomeString = 'incorrect';
      case 'dualrelease'
        outcomeString = 'dualrelease';
    end

    axH = subplot(spSz{:}, 1);						% default axes are 0 to 1
    
    set(axH, 'Visible', 'off');
    set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

    numTrials = nTrial;
    
    set(gcf, 'Visible', 'off'); % hide figure during text
                                % drawing - kludge

    text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 14);

    elMin = round((now - datenum(input.startDateVec)) * 24*60);
    startStr = datestr(input.startDateVec, 'HH:MM');
    text(0.00, 1.05, {'Subject:', 'Start time + elapsed:', 'Reward vol (m\pm std):'});
    text(0.60, 1.05, ...
         { sprintf('%2d', input.subjectNum), ...
           sprintf('%s + %2dm', ...
                   startStr, elMin), ...
           sprintf('%.1f s     \t(%g ms\\pm %g ms)', ...
                   nansum(juiceTimesSumMsV./1000), ...
                   chop(nanmean(juiceTimesSumMsV),2), ...
                   chop(nanstd(juiceTimesSumMsV),2)), ...
         });
    t2H(1) = text(0.00, 0.87, {'Trials:', 'Correct:', 'Early:', 'Failed:', 'Incorrect:', 'DualRelease:', 'MeanLeverLatency:'});
    t2H(2) = text(0.37, 0.87, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
                        sprintf('%d', nEarly), sprintf('%d', nIg), sprintf('%d', nInc), sprintf('%d', nDR)});
    t2H(3) = text(0.54, 0.87, {' ', sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
                        sprintf('%.0f%%', nEarly / numTrials * 100.0), ...
                        sprintf('%.0f%%', nIg / numTrials * ...
                                100.0)...
                        sprintf('%.0f%%', nInc / numTrials * 100.0), ...
                        sprintf('%.0f%%', nDR / numTrials * 100.0)});
    t2H(4) = text(0.88, 0.94, {'Left:', ...
                             sprintf('%2.0d', nLeftTrs), ...
                             sprintf('%2.0f%%', pLeftCorr*100 ), ...
                             sprintf('%2.0f%%', pLeftEarly*100), ...
                             sprintf('%2.0f%%', pLeftIg*100), ...
                             sprintf('%2.0f%%', pLeftIncorr*100), ...
                             sprintf('%2.0f%%', pLeftDR*100), ...
                             sprintf('%3.0f', leftLatencyMeanMs)});
    t2H(5) = text(1.16, 0.94, {'Right:', ...
                             sprintf('%2.0d', nRightTrs), ...
                             sprintf('%2.0f%%', pRightCorr*100 ), ...
                             sprintf('%2.0f%%', pRightEarly*100), ...
                             sprintf('%2.0f%%', pRightIg*100), ...
                             sprintf('%2.0f%%', pRightIncorr*100), ...
                             sprintf('%2.0f%%', pRightDR*100), ...
                             sprintf('%3.0f', rightLatencyMeanMs)});
    set(t2H, 'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left');

    
    % any funky block2 stuff here
    trPer80Str = regexprep(['[' num2str(input.trPer80V, '%3d') ']'], '[ *', '[');
    if input.rewardDoStaircase
        if input.rewardStaircaseDoSlidingWin
            sizeStr = sprintf('%d tr \\^ %g', ...
                              input.rewardStaircaseSlidingWinNTrials, ...
                              chop(input.rewardStaircaseAccelExponent,2));
        else % exponential weighting
            sizeStr = sprintf('%0.2f t weight', input.rewardStaircaseThisTrialWeight);
        end
        if input.rewardStaircaseDoAsymmetricTarget
            targStr = sprintf('%d', input.rewardStaircaseTargetUs/1000);
        else % symmetric
            targStr = '';
        end
        
        rewardStr = sprintf('Reward(stair min,t,max): %1d,%s,%1d ms; %s', ...
                            input.rewardStaircaseMinUs/1000, ...
                            targStr, ...
                            input.rewardStaircaseMaxUs/1000, ...
                            sizeStr);
    else
      rewardStr = sprintf('Reward(l,r) : %3.0f, %3.0f ms', ...
                    input.leftRewardUs/1000, ...
                    input.rightRewardUs/1000);
    end
    
    tStr = sprintf( ['Hold L,R; tf,dr: \t%3d, %3d; %3d, %3dms\n', ...
                     'Timeouts (e,m,i,dr):\t%4.1f, %4.1f, %4.1f, %4.1f s\n', ...
                     'React:\t%5.2f s;   ITI %5.0f ms\n', ...
                     '%s\n', ...
                     'SF: %3.2f cpd \n', ...                     
                     'trPer80: %s\n', ...
                     'Hold, React Median: %5.1f ms, %5.1f ms \n' ], ...
                    input.leftReqHoldMs, ...
                    input.rightReqHoldMs, ...
                    input.reactTooFastMs, ...
                    input.preventDualReleaseTimeMs, ...
                    input.timeoutEarlyMs/1000, ...
                    input.timeoutMissedMs/1000, ...
                    input.timeoutIncorrectMs/1000, ...
                    input.timeoutDualReleaseMs/1000, ...
                    input.reactTimeMs/1000, ...
                    input.itiFixedTimeMs, ...
                    rewardStr, ...
                    input.rightGratingSpatialFreqCPD, ...
                    trPer80Str, ...
                    median(holdV), ...
                    median(reactV(successIx)));

    text(0.0, 0.30, tStr, ...
         'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'top');
       
    set(gcf, 'Visible', 'on');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Hold time histogram
axH = subplot(spSz{:}, 4);
if input.reqHoldTimeMs{end} < 1000
    maxX = 2000;
else
    maxFail = input.reactTimeMs+ input.reqHoldTimeMs{end};
    maxX = ceil((maxFail+45)./500)*500;  % round up to nearest 500 ms.
end
if maxX > 10000, maxX = 10000; end

visIx = holdV <= maxX;
nVisPts = sum(visIx);
if nVisPts > 50
    binWidth = iqr(holdV(visIx)) ./ nVisPts.^(1/3);	% robust version of std
    nBins = ceil(maxX ./ binWidth);
else
    nBins = 10;
end

maxX=double(maxX);
edges = linspace(0, maxX, nBins);
Ns = histc(holdV(find(successIx)), edges);
Nf = histc(holdV(find(earlyIx)), edges);
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


if ~isempty(input.reactTooFastMs)
    yLim = get(gca, 'YLim');
    plot(double(input.reactTooFastMs) * [1 1], yLim, 'k--');
end
%%%%%%%%%%%%%%%%

%% 2 - Total React Time CDF
axH = subplot(spSz{:},5);

if ~isempty(leftReactTimes);
    cdfL= cdfplot(leftReactTimes(:));
    set(cdfL, 'Color', 'k');
    set(gca, 'XLim', [-1000 1000], ...
             'YLim', [0 1]);
    hold on;
end
if ~isempty(rightReactTimes);
    cdfR= cdfplot(rightReactTimes(:));
    set(cdfR, 'Color', 'b')
    set(gca, 'XLim', [-1000 1000], ...
             'YLim', [0 1]);
    hold on;
end
cdfTot=cdfplot(reactV);
set(cdfTot, 'Color', 'g');
set(gca, 'XLim', [-1000 1000], ...
         'YLim', [0 1]);
vH = vert_lines(0:200:1000);
set(vH, 'LineStyle', ':', 'Color', 0.5*[1 1 1]);
title('Total (g), Stim L (k), and R (b)');
xlabel('');
ylabel('F(x); RTs');


%% 3 - react time PDF
axH = subplot(spSz{:},11);
nPts = length(input.reactTimesMs);
visIx = reactV<=maxX;
nVisPts = sum(visIx);
if nVisPts > 50
    binWidth = iqr(reactV(visIx))./nVisPts.^(1/3); % robust version of std
    nBins = ceil(1000./binWidth);
else
    nBins = 10;
end
if nBins < 10; nBins=10; end

edges = linspace(-1000, 1000, nBins);
binSize = edges(2)-edges(1);

emptyIx = cellfun(@isempty, input.reactTimesMs);   % see above holdTimesM
if sum(emptyIx) > 0, input.reactTimesMs{emptyIx} = NaN; end
rV = reactV;

Ns = histc(rV(successIx), edges);
Nf = histc(rV(earlyIx), edges);
if sum(Ns)+sum(Nf) > 0
    bH = bar(edges+binSize/2, [Nf(:),Ns(:)], 'stacked');
    set(bH, 'BarWidth', 1, ...
            'LineStyle', 'none');
    cMap = get(gcf, 'Colormap');
    % flip colors, keep blue on top of red, note flipped in bar.m above
    set(bH(1), 'FaceColor', [0.6 0 0]);
    set(bH(2), 'FaceColor', [0 0 0.6]);      
end
title('Total React Times')
%% 4 - smoothed perf curve
axH = subplot(spSz{:},2);
hold on;
plot(smooth(double(successIx), ceil(nTrial/10), smoothType));

lH2 = plot(smooth(double(successIx), 100, smoothType));
set(lH2, 'Color', 'k', ...
         'LineWidth', 3);

lH3 = plot(smooth(double(ignoreIx), 100, smoothType));
set(lH3, 'Color', 'm');
lH4 = plot(smooth(double(earlyIx), 100, smoothType));
set(lH4, 'Color', 'c');
lH5 = plot(smooth(double(incorrectIx), 100, smoothType));
set(lH5, 'Color', 'g');
lH6 = plot(smooth(double(dualReleaseIx), 100, smoothType));
set(lH6, 'Color', [1 0.5 0.2]);

set([lH3 lH4 lH5 lH6], ...
    'LineWidth', 1, ...
    'LineStyle', '-');

anystack(lH6, 'bottom');
anystack(lH5, 'bottom');
anystack(lH4, 'bottom');
anystack(lH3, 'bottom');


ylabel('Percent Correct');
title('Total Perf Over Time')
set(gca, 'YLim', [0 1]);

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);


%% The left/right smoothed performance curve
axH = subplot(spSz{:},7);
hold on;
%plot(leftTrNs, smooth(double(leftCorr), ceil(nTrial/10), smoothType));
JH2 = plot(leftTrNs, smooth(double(leftCorr), 100, smoothType));
set(JH2, 'Color', 'k');
JH3 = plot(leftTrNs, smooth(double(leftIg), 100, smoothType));
set(JH3, 'Color', 'm');
JH5 = plot(leftTrNs, smooth(double(leftIncorr), 100, smoothType));
set(JH5, 'Color', 'g');
JH6 = plot(leftTrNs, smooth(double(leftDR), 100, smoothType));
set(JH6, 'Color', [1 0.5 0.2]);

set([JH2 JH3 JH5 JH6], ...
    'LineWidth', 2, ...
    'LineStyle', '-');

if ~isempty(leftTrNs);
    anystack(JH6, 'bottom');
    anystack(JH5, 'bottom');
    anystack(JH3, 'bottom');
end

hold on;
%plot(rightTrNs, smooth(double(rightCorr), ceil(nTrial/10), smoothType));
KH2 = plot(rightTrNs, smooth(double(rightCorr), 100, smoothType));
set(KH2, 'Color', 'k');
KH3 = plot(rightTrNs, smooth(double(rightIg), 100, smoothType));
set(KH3, 'Color', 'm')
KH5 = plot(rightTrNs, smooth(double(rightIncorr), 100, smoothType));
set(KH5, 'Color', 'g');
KH6 = plot(rightTrNs, smooth(double(rightDR), 100, smoothType));
set(KH6, 'Color', colOrange);

set([KH2 KH3 KH5 KH6], ...
    'LineWidth', 2, ...
    'LineStyle', '-.');

if ~isempty(rightTrNs);
    anystack(KH6, 'bottom');
    anystack(KH5, 'bottom');
    anystack(KH3, 'bottom');
end

% pct corr for seen stim
smPm = { 50, 'lowess' };
rC = smooth(double(rightCorr), smPm{:});
rI = smooth(double(rightIncorr), smPm{:});
lC = smooth(double(leftCorr), smPm{:});
lI = smooth(double(leftIncorr), smPm{:});

rCorrOnStim = rC ./ (rC+rI);
lCorrOnStim = lC ./ (lC+lI);
scH1 = plot(rightTrNs, smooth(rCorrOnStim, 50, 'lowess'), 'b-.');
scH2 = plot(leftTrNs, smooth(lCorrOnStim, 50, 'lowess'), 'b');

ylabel('Percent Correct');
set(gca, 'YLim', [0 1]);
title('L (solid) and R (dotted) Performance')

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);

%% Hold Start Difference graph

axH = subplot(spSz{:},6);
hold on;

hSDiffsRealSec = double(diff(holdStarts)/1000);
xs = 1:length(hSDiffsRealSec);
pH1 = plot(xs, cumsum(hSDiffsRealSec)./60, '.-');
maxMin = sum(hSDiffsRealSec)./60;

title('Time Working (min) vs. Trials');
%set(gca, 'DataAspectRatio', [100 10 1]);
xLim = trXLim;
set(gca, 'XLim', xLim);
% compute data aspect ratio manually: matlab will keep y dim of
% plot box fixed and change x and we want the reverse
yLim = [0 max(xLim(2)/10, maxMin*1.1)];   
set(gca, 'YLim', yLim);
pH2 = plot(xLim, xLim/10, 'k--');

%% The switch/%left percent plot
axH = subplot(spSz{:},3);
hold on;

trNs = 1:nTrial;
leftRew = juiceTimesSumMsV(leftIx & successIx);
lTrNs = trNs(leftIx & successIx);
rightRew = juiceTimesSumMsV(~leftIx & successIx);
rTrNs = trNs(~leftIx & successIx);

% setup plot
[axH cH1 cH2] = plotyy(0,1,0,1);

% compute percent left in a window like the staircase window
if input.rewardDoStaircase
    if input.rewardStaircaseDoSlidingWin
        swY = smooth(tLeftTrial(~earlyIx), 'moving', input.rewardStaircaseSlidingWinNTrials*2);
    else
        swY = smooth(tLeftTrial(~earlyIx), 15, smoothType);
    end
else
    swY = smooth(tLeftTrial(~earlyIx), 50, smoothType);
end

% plot trial fractions
axis(axH(1));
hold on;
zH = plot(trNs(~earlyIx), swY);
zSwH = plot(trNs, smooth(tSwitchV, 50, smoothType));
set(zH, 'Color', colYellow);
set(zSwH, 'Color', colPurple, ...
          'LineStyle', ':', ...
          'LineWidth', 2);
%set([zH zSwH], ...
%    'LineWidth', 2);


% find resp dir and plot
respS = tFirstReactReleaseIsLeft;  % start with any trial rel L or R first
respS(missIx | dualReleaseIx | earlyIx) = NaN;  % blank out miss/dr/early, keeping corr, incorr

negTRIx=(respS==-1);
respS(negTRIx)= NaN;

resH= plot(trNs, smooth(respS, 20, smoothType));
set(resH, 'Color', colOrange, ...
              'LineWidth', 1);

% find resp switch and plot
switchRespIx = ~isnan(respS);
switchVals = [NaN diff(respS(switchRespIx)) ~= 0];
switchS = NaN*respS;
switchS(switchRespIx) = switchVals;

swH= plot(trNs, smooth(switchS, 50, smoothType));
set(swH, 'Color', colPurple, ...
              'LineWidth', 2);

ylabel('% Left/Switch');
anystack([resH, swH], 'bottom');
anystack([zH, zSwH], 'bottom')

% plot reward amount
set(axH(2), 'NextPlot', 'add');
if ~isempty(lTrNs)
    rewLH = plot(axH(2), lTrNs, leftRew);
else  
    rewLH = [];
end
if ~isempty(rTrNs);
    rewRH = plot(axH(2), rTrNs, rightRew);
else
    rewRH = [];
end
set([rewLH rewRH], 'LineStyle', 'none', ...
                  'Marker', '.');
set(rewLH, 'Color', 'b');
set(rewRH, 'Color', 'k');

set(axH(1), 'YColor', 'k');
set(axH(2), 'YColor', 'b');

if input.rewardDoStaircase
    yT = [0 ...
          input.rewardStaircaseMinUs ...
          input.rewardStaircaseTargetUs ...
          input.rewardStaircaseMaxUs ] / 1000;
else
    yT = [0 ...
          input.leftRewardUs ...
          input.rightRewardUs] / 1000;
end
yT = unique(yT);
yTL = sprintf_vector('%d', yT);
if ~input.rewardStaircaseDoAsymmetricTarget
    yTL{3} = '';
end

y1T = [0 1];
y1T = sort([y1T input.stimProbAvgLeft]);

yL = get(axH(2), 'YLim');
yL(1) = 0; 
yL(2) = nanmax([leftRew(:); rightRew(:); yT(end)*1.03]);
set(axH(2), ...
  'YLim', yL, ...
  'YTickLabel', yTL, ...
  'YTick', yT);
biasFoundIx = ~isnan(tRewardRunningLeftBias);
biasH = plot(trNs(biasFoundIx), tRewardRunningLeftBias(biasFoundIx));
set(biasH, 'Color', colOrange, ...
  'LineWidth', 2);

set(axH, 'NextPlot', 'add');
set(axH(1), 'YLim', [0 1], ...
            'YTick', y1T);
title('Stim: L(y) Sw(p:) Resp: L(o) lBias(O) Sw(p)')
set(axH, 'XLim', [0 nTrial]);
set(axH(2), 'Position', get(axH(1), 'Position'));
get(axH, 'XLim');
get(axH, 'XLimMode');
hold(axH(2), 'off')

%%%%%%%%%%%%%%%%
%% subplot: time between trials (sp 8)

axH = subplot(spSz{:},8);
hold on;
holdStarts = double(cellvect2mat_padded(input.holdStartsMs));
leftHoldStarts = holdStarts(leftIx);
rightHoldStarts = holdStarts(~leftIx);
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

%% Left/Right Mean Hold/React plot
axH = subplot(spSz{:}, 10);
hold on;

vL1 = smooth(leftHolds, 50, 'rloess');
vR1 = smooth(rightHolds, 50, 'rloess');
vyL1 = smooth(leftReacts, 35, 'rloess');
vyR1 = smooth(rightReacts, 35, 'rloess');
[axH pH1 pH2] = plotyy(1,1,1,1);
set(axH, 'NextPlot', 'add');

% first axes
if ~isempty(vL1);
    hH(1) = plot(axH(1), leftECTrs, vL1, 'LineWidth', 1.5);
    
    c1 = 'k';
    set(hH, 'Color', c1);
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
    ylabel('Hold Times (ms) - Corr+Early');
end
if ~isempty(vR1);
    hH(1) = plot(axH(1), rightECTrs, vR1, '-.', 'LineWidth', 1.5);
    
    c1 = 'k';
    set(hH, 'Color', c1);
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
    ylabel('Hold Times (ms) - resp+early');
end

% 2nd axes
if ~isempty(vyL1);
    hyH(1) = plot(axH(2), leftCTrs, vyL1, 'LineWidth', 1.5);    
    c2 = 'b';
    set(hyH, 'Color', c2);    
    %    set(axH(2), 'YLimMode', 'auto', ...
    %                'YTickMode', 'auto', ...
    %                'YTickLabelMode', 'auto', ...
    %                'YColor', c2)
    ylabel(axH(2), 'React time (ms) - corr');
end
if ~isempty(vyR1);
    hyH(1) = plot(axH(2), rightCTrs, vyR1, '-.', 'LineWidth', 1.5);    
    c2 = 'b';
    set(hyH, 'Color', c2);    
    set(axH(2), 'YLimMode', 'auto', ...
                'YTickMode', 'auto', ...
                'YTickLabelMode', 'auto', ...
                'YColor', c2)
    ylabel(axH(2), 'React time (ms) - resp made');
end
title('L (solid) and R (dotted) resp');

%% Left/Right Latency Plot
axH = subplot(spSz{:}, 12);
hold on;

desLTrIx = (successIx | incorrectIx | dualReleaseIx) & tLeftTrial==1;
desRTrIx = (successIx | incorrectIx | dualReleaseIx) & tLeftTrial==0;


trNs = 1:nTrial;
latencyV = leverLatencyRPosMs;
latencyV(latencyV > 250) = 250;
latencyV(latencyV < -250) = -250;
latencyV(isnan(latencyV) & tFirstReactReleaseIsLeft==1) = -250;
latencyV(isnan(latencyV) & tFirstReactReleaseIsLeft==0) = +250;

leftStimLatencyMs = latencyV(desLTrIx);
rightStimLatencyMs = latencyV(desRTrIx);

pH3=plot( trNs(desLTrIx&dualReleaseIx), latencyV(desLTrIx & dualReleaseIx), 'b+');
pH4=plot( trNs(desRTrIx&dualReleaseIx), latencyV(desRTrIx & dualReleaseIx), 'k+');

pH1 = plot( trNs(desLTrIx), smooth(leftStimLatencyMs, 25, 'loess'), 'b');
pH2 = plot( trNs(desRTrIx), smooth(rightStimLatencyMs, 25, 'loess'), 'k');

pH5=plot( trNs(desLTrIx & ~dualReleaseIx), latencyV(desLTrIx & ~dualReleaseIx), 'b.');
pH6=plot( trNs(desRTrIx & ~dualReleaseIx), latencyV(desRTrIx & ~dualReleaseIx), 'k.');

%symmetric axis
%yLim = get(gca, 'YLim');
%yLim = max(abs(yLim))*[-1 1];
yLim = [-260 260];
set(gca, 'YLim', yLim);
% dashed zero line
xLim = [0 nTrial];
lH = plot(xLim, [0 0], 'r--');
anystack(lH, 'bottom');
set(gca, 'XLim', xLim, ... 
         'YDir', 'reverse');

set([pH1 pH2], 'LineWidth', 2);
set([pH3; pH4], 'MarkerSize', 4); 
box on;
%set(gca, 'YGrid', 'on');

% annotate
ylabel('mean lever latency (R +)');
title('latency on stim: R + k; L - b');

%% add marginal distribution to right
lastPos = get(axH, 'Position');
axH = axes('Units', 'Normalized',...
           'Position', [0.91 lastPos(2) 0.10 lastPos(4)]); 
hold on;

edges = [-255:10:265];


NL = histc(leftStimLatencyMs, edges);
NR = histc(rightStimLatencyMs, edges);
if isempty(NL)
    NL = edges*0;
end
if isempty(NR)
    NR = edges*0;
end


%desNs = 2:(length(edges)-1);
desNs = 1:length(edges);
lV = NL(desNs) ./ max(NL(desNs));
rV = NR(desNs) ./ max(NR(desNs));
endNs = [1 length(edges)];
capVal = 0.4;
capLIx = lV(endNs) >= capVal;
capRIx = rV(endNs) >= capVal;
if any(capLIx), lV(endNs(capLIx)) = capVal; end
if any(capRIx), lV(endNs(capRIx)) = capVal; end

tX = edges(desNs);
lH1=plot([0 0], [min(tX) max(tX)], '-k');  % 0 line for reference
set(lH1, 'Color', 'r');

lV = smooth(lV, 5, 'lowess');
rV = smooth(rV, 5, 'lowess');
lH1 = plot(lV, tX, 'b');
lH2 = plot(rV, tX, 'k');

if any(capLIx), plot(lV(capLIx), edges(capLIx), 'b+'); end
if any(capRIx), plot(rV(capRIx), edges(capRIx), 'k+'); end

set([lH1 lH2], 'LineWidth', 2);

set(gca, 'XLim', [-.05 1], ...
         'YLim', yLim, ... % from above
         'Visible', 'off', ...
         'YDir', 'reverse');


%% Add a save button
outDir = cs.behavPdfPath;
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
        
        trPer80Changed = false;
        for iR = 1:nRows
            tVarName = tChangedList{iR,1};
            tChangedTo = tChangedList{iR,2};
            % special case odds changes
            if strfind(tVarName, 'trPer80Level')
                trPer80Changed = true;
            else
                changedStrOut{end+1} = sprintf('Tr %3d - %s: -> %g', ...
                                               tDesN, tVarName, tChangedTo);
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
end

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
