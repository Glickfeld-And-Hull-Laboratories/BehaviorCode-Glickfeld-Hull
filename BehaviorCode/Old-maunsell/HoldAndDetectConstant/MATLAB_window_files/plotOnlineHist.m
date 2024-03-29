function input = plotOnlineHist(data_struct, input)

% essential consts
figNum = 4;
name = 'HoldAndDetectConstant';
cs = holdanddetect_constants;
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
  figPos = [700 140 780 770];
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
if isfield(input, 'tDoNoStimulusChange'),
    dnscIx = celleqel2mat_padded(input.tDoNoStimulusChange)==1;
    dnscTrs = find(dnscIx);
    dnscCorrIx = dnscIx & successIx;
    dnscCorrTrs = find(dnscCorrIx);
    dnscCorrIx = dnscCorrIx(dnscIx);
end

%  TODO:  this should be automated in a loop, or better should not have to convert from cell at all
tGratingDurationMsV = celleqel2mat_padded(input.tGratingDurationMs);
tGratingDirectionDeg = celleqel2mat_padded(input.tGratingDirectionDeg);
tGratingContrast = celleqel2mat_padded(input.tGratingContrast, NaN, 'double');
tBaseGratingContrast = celleqel2mat_padded(input.tBaseGratingContrast, NaN, 'double');
tGratingSpeedDPS = celleqel2mat_padded(input.tGratingSpeedDPS, NaN, 'double');
tBaseGratingSpeedDPS = celleqel2mat_padded(input.tBaseGratingSpeedDPS, NaN, 'double');
tLaserPowerMw = celleqel2mat_padded(input.tLaserPowerMw);
tLaserRampLengthMsV = celleqel2mat_padded(input.tLaserRampLengthMs);
tLaserDoLinearRampV = celleqel2mat_padded(input.tLaserDoLinearRamp);
tLaserDoPulseTrainV = celleqel2mat_padded(input.tLaserDoPulseTrain);
tTrialLaserPowerMwV = celleqel2mat_padded(input.tTrialLaserPowerMw);
tBlock2TrialNumberV = celleqel2mat_padded(input.tBlock2TrialNumber);
tStimulusNumberV = celleqel2mat_padded(input.tStimulusNumber);
tLaserBaselinePowerMwV = celleqel2mat_padded(input.tLaserBaselinePowerMw);

reactV = celleqel2mat_padded(input.reactTimesMs);
holdV = celleqel2mat_padded(input.holdTimesMs);  % sometimes on client restart have empty elements here; should figure out why
juiceTimesMsV = cellfun(@sum, input.juiceTimesMsCell);
juiceTimesMsV(juiceTimesMsV==0) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% figure out the number of stimuli
lPowerV = tLaserPowerMw;
if all(lPowerV) == 0
  nLP = 0;
else
  nLP = length(chop(unique(lPowerV(~isnan(lPowerV))),4));
end
if input.doContrastDetect
    vPowerV = abs(tGratingContrast-tBaseGratingContrast)*100;
elseif input.doOriDetect
    vPowerV = tGratingDirectionDeg;
elseif input.doSpeedDetect
    vPowerV = abs(tGratingSpeedDPS-tBaseGratingSpeedDPS);
end

if all(vPowerV) == 0
  nVP = 0;
else
  nVP = length(chop(unique(vPowerV(~isnan(vPowerV))),4));
end

if nVP >= 1
  showLaserStim = false;
  nStims = nVP;
elseif nLP >= 1
  showLaserStim = true;
  nStims = nLP;
else
  if all(isnan(vPowerV))
    showLaserStim = true;
    nStims = 1;
  elseif all(isnan(lPowerV))
    showLaserStim = false;
    nStims = 1;
  elseif all(~isnan(vPowerV) & ~isnan(lPowerV))
    % both laser and stim each trial
    showLaserStim = false;
    nStims = 1;
  else
    disp(sprintf('nLaser %d, nVis %d', nLP, nVP));
    error('confused on number of laser and vis levels: contrast/dir base and step set?');
  end
end

%% disp
%stimFrames = cellfun(@length,input.announceStimulusTimes.driftStimulus)

% misc error checking
%assert(~input.doLaserStim || all(tLaserDoLinearRampV | tLaserDoPulseTrainV));

% figure out block2 stim levels
showBlock2 = input.doBlock2;
if showBlock2
  %assert(all(tBlock2TrialNumberV == 0 | tBlock2TrialNumberV == 1), 'xml bug?');

  if input.block2DoGratingAppearance
    block2Name = 'Grating appearance';
    block2LevelNames = sprintf_vector('grating %d', [1 2]);  % should improve this

    block2Tr1Ix = tBlock2TrialNumberV == 0;
    block2Tr2Ix = tBlock2TrialNumberV == 1;
  
  elseif input.block2DoRampLength
    block2Levels = [input.laserRampLengthMs input.block2RampLengthMs2];
    block2Name = 'Ramp length';
    block2LevelNames('%dms', block2Levels);
    b2Temp = unique(tLaserRampLengthMsV);
    assert(length(b2Temp<2) || all(sort(block2Levels) == b2Temp), 'must reset after changing block2 levels');
    nL = length(block2Levels);
    assert(nL <= 2);
  
    block2Tr1Ix = tLaserRampLengthMsV == block2Levels(1);
    block2Tr2Ix = tLaserRampLengthMsV == block2Levels(2);
  elseif input.block2DoRampVTrain
    block2Levels = [1 0];
    block2Name = 'Ramp v. train';
    block2LevelNames = {'ramp', 'train'};
    nL = 2;
    
    block2Tr1Ix = tLaserDoLinearRampV == 1;
    block2Tr2Ix = tLaserDoLinearRampV == 0;
  elseif input.block2DoTrialLaser
    block2Levels = [0 1];
    block2Name = 'trial laser';
    block2LevelNames = {'b1', 'b2'};
    nL = 2;
    
    block2Tr1Ix = tBlock2TrialNumberV == 0;
    block2Tr2Ix = tBlock2TrialNumberV == 1;
  else
    error('unable to find block2 variable!');
  end
end

% figure out misc trial types
isRampTrial =  input.tLaserDoLinearRamp{input.trialSinceReset};
isTrainTrial =  input.tLaserDoPulseTrain{input.trialSinceReset};
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

        axH = subplot(spSz{:}, 1);            % default axes are 0 to 1
        
  set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

        numTrials = nTrial;
        elMin = round((now - datenum(input.startDateVec)) * 24*60);
        startStr = datestr(input.startDateVec, 'HH:MM');
        
        set(gcf, 'Visible', 'off'); % hide figure during text
                                    % drawing - kludge

        text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 14);
        text(0.00, 1.05, {'Subject:', 'Start time + elapsed:', 'Reward vol (m \pm std):'});
        text(0.60, 1.05, ...
             { sprintf('%2d', input.subjectNum), ...
               sprintf('%s + %2dm', startStr, elMin), ...
               sprintf('%.1fs     (%gms \\pm %gms)', ...
                       nansum(juiceTimesMsV./1000), ...
                       chop(nanmean(juiceTimesMsV),2), ...
                       chop(nanstd(juiceTimesMsV),2))
             });
  t2H(1) = text(0.00, 0.9, {'Trials:', 'Correct:', 'Early:', 'Failed:'});
  t2H(2) = text(0.35, 0.9, {sprintf('%d', nTrial), sprintf('%d', nCorr), ...
                            sprintf('%d', nFail), sprintf('%d', nIg)});
  t2H(3) = text(0.54, 0.9, {' ', ...
        sprintf('%.0f%%', nCorr / numTrials * 100.0), ...
        sprintf('%.0f%%', nFail / numTrials * 100.0), ...
        sprintf('%.0f%%', nIg / numTrials * 100.0)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');
        tStr = sprintf(['Hold, react median:\n', ...
                        '   %5.1f ms, %5.1f ms\n', ...
                        '\n'], ...
                       median(holdV), ...
                       median(reactV(successIx)));
               
        text(0.8, 0.9, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left');

        if input.laserRampDoExpRamp == 0, 
          rampShapeStr = 'lin'; 
        else
          rampShapeStr = 'exp';
        end
        
        trPer80Str = regexprep(['[' num2str(input.trPer80V, '%3d') ']'], '[ *', '[');
        if input.doBlock2SeparateOdds
            block2TrPer80Str = regexprep(['[' num2str(input.block2TrPer80V, '%3d') ']'], '[ *', '[');
        else
            block2TrPer80Str = '';
        end

        % any funky block2 stuff here
        if showBlock2
          if input.block2DoRampLength
            block2Str = sprintf('ramp length: %d %d ms', input.laserRampLengthMs, input.block2RampLengthMs2);
          elseif input.block2DoGratingAppearance
            block2Str = sprintf('grating 1: %s\n          grating 2: %s', subPpGrating(input), subPpGrating(input,true));
          elseif input.block2DoRampVTrain
            block2Str = sprintf('ramp v train: ramp %d (%s)+%d ms, bl %g mW\n  train, pulse %d period %d dur %d, bl %g mW', ...
                                input.laserRampLengthMs, ...
                                rampShapeStr, ...
                                input.laserRampExtraConstantLengthMs, ...
                                chop(input.block2RvtRampBaselinePowerMw,2), ...
                                input.laserPulseLengthMs, input.laserPulsePeriodMs, ...
                                input.laserTrainLengthMs, ...
                                chop(input.block2RvtTrainBaselinePowerMw,2));
          elseif input.block2DoTrialLaser
            block2Str = sprintf('trial laser 1: %s\n          trial laser 2: %s', subPpTrialLaser(input,1), subPpTrialLaser(input,2));
          end
        else
          block2Str = 'off';
        end

        % stimulus str
        if input.doBlock2
          if input.doVisualStim & input.block2DoTrialLaser             
            stimStr = subPpGrating(input);
          else
            stimStr = 'block2 controlled';
          end
        else
          stimStr = '';
          if input.doLaserStim
            stimStr = subPpLaser(input);
          end
          if input.doVisualStim
            stimStr = strcat(stimStr, subPpGrating(input));
          end
          if isfield(input, 'doTrialLaser') && input.doTrialLaser
            stimStr = strcat(stimStr, subPpTrialLaser(input,1));
          end
        end
        if ~input.doBlock2
          if input.laserBaselinePowerMw > 0
            stimStr = strcat(stimStr, sprintf(', baseline %gmW', chop(input.laserBaselinePowerMw,2)));
          end
        end
        
        if input.doGeomHoldDist
          randTypeStr = 'geo';
          geoExtraStr = sprintf('(geo m%d)', input.geomHoldMeanMs);
        else 
          randTypeStr = 'unif';
          geoExtraStr = '';
        end

        if input.doBlock2 && input.doBlock2SeparateReward
            b2RewardStr = sprintf('(b2: %d,%d ms)', ...
                                  round(input.block2MinRewardUs/1000), ...
                                  round(input.block2MaxRewardUs/1000));
        else
            b2RewardStr = [];
        end
            
        tStr = sprintf( ['Hold (f+r/%s, tf):       %3d+%3d,%3dms %s\n', ...
                         'Timeouts (e,m):        %4.1f, %4.1fs\n', ...
                         'React: %5.2fs;           ITI: %dms\n', ...
                         'Reward min, max:      %3.0f, %3.0fms   %s\n', ...
                         '\n', ...
                         'trPer80: %s\n', ...
                         'block2TrPer80: %s\n', ...
                         'Stim: %s\n', ...
                         'block2: %s\n' ], ...
                        randTypeStr, ...
                        input.fixedReqHoldTimeMs, ...
                        input.randReqHoldMaxMs, ...
                        input.tooFastTimeMs, ...
                        geoExtraStr, ...
                        input.earlyTimeoutMs/1000, ...
                        input.missedTimeoutMs/1000, ...
                        input.reactTimeMs/1000, ...
                        input.itiTimeMs, ...
                        input.minRewardUs/1000, ...
                        input.maxRewardUs/1000, ...
                        b2RewardStr, ...
                        trPer80Str, ... 
                        block2TrPer80Str, ...
                        stimStr, ...
                        block2Str);

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

axH = subplot(spSz{:}, 3);
if input.fixedReqHoldTimeMs+input.randReqHoldMaxMs < 1000
  maxX = 2000;
else
  maxFail = double(input.reactTimeMs+input.fixedReqHoldTimeMs+input.randReqHoldMaxMs);
  maxX = ceil((maxFail+45)./500)*500;  % round up to nearest 500 ms.
end

visIx = holdV <= maxX;
nVisPts = sum(visIx);
if nVisPts > 50
  binWidth = iqr(holdV(visIx)) ./ nVisPts.^(1/3); % robust version of std
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

title('total hold times');


if ~isempty(input.tooFastTimeMs)
  yLim = get(gca, 'YLim');
  plot(double(input.tooFastTimeMs) * [1 1], yLim, 'k--');
end
%%%%%%%%%%%%%%%%

%% 2 - react time CDF
axH = subplot(spSz{:},4);
reacts = [input.reactTimesMs{:}];
cdfplot(reacts);
set(gca, 'XLim', [-1000 1000], ...
         'YLim', [0 1]);
hold on;
if isfield(input, 'tDoNoStimulusChange'),
    if sum(dnscIx)>0,
        dnscP = cdfplot(reacts(dnscIx));
        set(dnscP,'Color', [0 0.6 0]);
        set(gca, 'XLim', [-1000 1000], 'YLim', [0 1]);
    end
end
vH = vert_lines(0:200:1000);
set(vH, 'LineStyle', ':', 'Color', 0.5*[1 1 1]);
title('');
xlabel('');

%%%%%%%%%%%%%%%%

%% 3 - react time PDF
axH = subplot(spSz{:},7);
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
title('reaction times');
%%%%%%%%%%%%%%%%

%% 4 - smoothed perf curve
axH = subplot(spSz{:},2);
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
     
if isfield(input, 'tDoNoStimulusChange') & sum(dnscIx)>0,
    lH5 = plot(dnscTrs, smooth(dnscCorrIx, 100, smoothType));
    set(lH5, 'Color',[0 0.6 0], ...
           'LineWidth', 2);
end
if showBlock2
  b2T1Ns = find(block2Tr1Ix);
  b2T2Ns = find(block2Tr2Ix);
  lH4a = plot(b2T1Ns, smooth(double(failureIx(block2Tr1Ix)), 75, smoothType));
  lH4b = plot(b2T2Ns, smooth(double(failureIx(block2Tr2Ix)), 75, smoothType));
  set([lH4a lH4b], 'LineWidth', 2, 'LineStyle', '-.');
  set(lH4a, 'Color', 0.8*[0 1 1]);
  set(lH4b, 'Color', 0.8*[1 1 0]);
  anystack([lH4a lH4b], 'bottom');
else
  lH4 = plot(smooth(double(failureIx), 100, smoothType));
  set(lH4, 'Color', 0.8*[0 1 1], ...
           'LineWidth', 2, ...
           'LineStyle', '-.');
  anystack(lH4, 'bottom');
end

anystack(lH3, 'bottom');


ylabel('pct correct');
set(gca, 'YLim', [0 1]);

trXLim = [0 nTrials]; %get(gca, 'XLim');
set(gca, 'XLim', trXLim);
    
      

%%%%%%%%%%%%%%%%

%% 6 - trial length plot
axH = subplot(spSz{:},6);
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
set(axesH, 'XLim', xLim);
lH = plot(xLim, 20*[1 1], '--k');

set(axesH(1),'YLim', [0 121], ...
             'YTick', [0 40 80 120], ...
             'YTickLabelMode', 'auto', ...
             'YColor', 'k');

if ~isempty(axesH(2)) 
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
axH = subplot(spSz{:},5);
hold on;

hSDiffsRealSec = diff(holdStarts)/1000;
xs = 1:length(hSDiffsRealSec);
pH1 = plot(xs, cumsum(hSDiffsRealSec)./60, '.-');
maxMin = sum(hSDiffsRealSec)./60;

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
    if isfield(input, 'tDoNoStimulusChange') && sum(dnscIx)>0,
        dnscCorrIx = dnscIx & successIx;
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
    ylabel(axH(2), 'React time (ms) - corr');
end

title('mean react/hold over time (robust fit)');

% used to have % 9 - reward sizes over correct trials


%% 10 -  number of responses during react window - only for 1 stim
%level
if nStims == 1
  axH = subplot(spSz{:}, 10);
    winLen = 400;
  if median(reactV(successIx)) <= 400;
    winStart = 150;
  else
    winStart = 250;
  end
  winReal = reactV > winStart & reactV <= (winStart+winLen);
  baselineRange = [-winLen*1.5 winLen+winStart];         % From -600ms --> 550ms/650ms
  nBaseWins = range(baselineRange)./winLen;
  winBaseline = reactV > baselineRange(1) & reactV <= baselineRange(2);
  winBefore = winBaseline ./ nBaseWins;
  if nTrial <= 80
    reactionRate = (sum(winReal)./nTrial)*100;
  else
    recentTrials = winReal(nTrial-79:nTrial);                         % Last 80 Trials
    reactionRate = (sum(recentTrials)./length(recentTrials))*100;
  end

  %% Calculates % of responses within a specific range out of total trials - Blue 
  %% and trials near stimulus change, from 600ms before to end of react window - Black
  
  hold on;
  p1H = plot(smooth(double(winReal), ceil(nTrial/10), 'lowess'));
  set(p1H, 'LineWidth', 3, 'Color', [0.25 0.3 0.9]);
  nSm = 5;
  v1 = smooth(double(winReal), ceil(nTrial./nSm), 'lowess');
  v2 = smooth(double(winBefore), ceil(nTrial./nSm), 'lowess');
  v0 = v1 ./ (v1+v2);
  p2H = plot(v0, 'k');
  set(p2H, 'LineWidth', 1, 'LineStyle', '--');
  %p2H = plot(smooth(double(winBefore), ceil(nTrial/10), 'lowess'));
  %p3H = plot(smooth(double(winAfter), ceil(nTrial/10), 'lowess'));

  title(['Reaction to Stimulus (last 80): ',num2str(reactionRate,3),'%'])
  if winStart == 150;
    ylabel('% Correct in 150-550ms')
  elseif winStart == 250;
    ylabel('% Correct in 250-650ms')
  end
  %xLim = [0 nTrials];
  xLim = trXLim;
  set(gca, 'YLim', [0 1], ...
           'XLim', xLim);
  plot(xLim, 0.30*[1 1], '--','Color', [0.65 0.55 0.95])
  plot(xLim, 0.55*[1 1], '--','Color',[0.78 0.72 0.00]);
  
  axH = subplot(spSz{:}, 4);
  hold on;
  v2H = vert_lines([winStart winStart+winLen]);
  set(v2H, 'Color', 'g', 'LineStyle', '--');
end

%%%%%%%%%%%%%%%%

%% 11 - plot psych curves  (do for all)
if nStims >= 1
  axH = subplot(spSz{:}, 11);
  hold on;

  if showLaserStim
    stimPowerV = chop(lPowerV,2);
  else
    stimPowerV = chop(vPowerV,2);
  end
  powerLevels = unique(stimPowerV(~isnan(stimPowerV)));
  nL = length(powerLevels);
  
  % init structure for later use in for loop
  pFNames = {'correct','early','total', ...
             'correct1','early1','total1', ...
             'correct2','early2','total2', ...
             'correctRecent','earlyRecent','totalRecent'};
  paramC = cell(1,length(pFNames)*2);
  paramC(1:2:end) = pFNames(:);
  [paramC{2:2:end}] = deal(NaN*ones([1 nL]));
  powerP = struct(paramC{:});

  maxR = max(input.fixedReqHoldTimeMs)+max(input.randReqHoldMaxMs);
  minR = input.fixedReqHoldTimeMs;
  halfR = (maxR-minR)/2 + minR;
  firstHalfIx = reqHoldV <= halfR;
  secondHalfIx = reqHoldV > halfR;
  lastFewTrIx = false(size(firstHalfIx));  
  lastFewTrIx(find(successIx | missIx, 50, 'last')) = true;


  for iP=1:nL
    tP = powerLevels(iP);
    trIx = stimPowerV == tP;
    
    powerP(1).correct(iP) = sum(trIx & successIx);
    powerP.early(iP) = sum(trIx & earlyIx);
    powerP.total(iP) = sum(trIx);

    powerP.correct1(iP) = sum(trIx & successIx & firstHalfIx);
    powerP.early1(iP) = sum(trIx & earlyIx & firstHalfIx);
    powerP.total1(iP) = sum(trIx & firstHalfIx);

    powerP.correct2(iP) = sum(trIx & successIx & secondHalfIx);
    powerP.early2(iP) = sum(trIx & earlyIx & secondHalfIx);
    powerP.total2(iP) = sum(trIx & secondHalfIx);
  
    powerP.correctRecent(iP) = sum(trIx & successIx & lastFewTrIx);
    powerP.earlyRecent(iP) = sum(trIx & earlyIx & lastFewTrIx);
    powerP.totalRecent(iP) = sum(trIx & lastFewTrIx);

    if showBlock2
      powerP.block2Correct1(iP) = sum(trIx & block2Tr1Ix & successIx);
      powerP.block2Correct2(iP) = sum(trIx & block2Tr2Ix & successIx);
      powerP.block2Early1(iP) = sum(trIx & block2Tr1Ix & earlyIx);
      powerP.block2Early2(iP) = sum(trIx & block2Tr2Ix & earlyIx);
      powerP.block2Total1(iP) = sum(trIx & block2Tr1Ix);
      powerP.block2Total2(iP) = sum(trIx & block2Tr2Ix);
    end
    
  end
  powerP.powerRecentIx = ismember(powerLevels, unique(stimPowerV(lastFewTrIx)));

  pctCorr = (powerP.correct) ./ (powerP.total-powerP.early) ...
            *100;
  
  pctCorr1 = (powerP.correct1) ./ (powerP.total1-powerP.early1) * 100;
  pctCorr2 = (powerP.correct2) ./ (powerP.total2-powerP.early2) * 100;
  pctCorrRecent = (powerP.correctRecent) ./ (powerP.totalRecent-powerP.earlyRecent) * 100;
  if showBlock2
    pctCorr1Block2 = powerP.block2Correct1 ./ (powerP.block2Total1-powerP.block2Early1) * 100;
    pctCorr2Block2 = powerP.block2Correct2 ./ (powerP.block2Total2-powerP.block2Early2) * 100;
    des1Ix = ~isnan(pctCorr1Block2);
    des2Ix = ~isnan(pctCorr2Block2);
    pH=[NaN NaN];
    if sum(des1Ix)>0
      pH(1) = plot(powerLevels(des1Ix), pctCorr1Block2(des1Ix), '.-b');
    end
    if sum(des2Ix)>0
      pH(2) = plot(powerLevels(des2Ix), pctCorr2Block2(des2Ix), '.-');
      set(pH(2), 'Color', yColor);
    end
  else
    pH = plot(powerLevels, pctCorr, '.-b');
  end
  if all(~isempty(pH))
    set(pH, 'LineWidth', 2)
  end

  pH1 = plot(powerLevels, pctCorr1, '-g');
  pH2 = plot(powerLevels, pctCorr2, '-r');
  if showBlock2
    set([pH1 pH2], 'LineStyle', 'none', ...
                   'Marker', '+');
  end

  if showBlock2
    pH3 = plot(powerLevels(powerP.powerRecentIx), pctCorrRecent(powerP.powerRecentIx), 'xk');
  else
    pH3 = plot(powerLevels(powerP.powerRecentIx), pctCorrRecent(powerP.powerRecentIx), '-k');
  end
  
  % baseline annotate
  bH1 = []; bH2 = []; bH = [];
  if showBlock2
    unqBase1 = unique(tLaserBaselinePowerMwV(block2Tr1Ix));
    unqBase2 = unique(tLaserBaselinePowerMwV(block2Tr2Ix));
    if any(unqBase1>0)
      bH1 = plot(unqBase1, 5+unqBase1*0);
      set(bH1, 'Color', 'b', 'MarkerFaceColor', 'b');
    end
    if any(unqBase2>0)
      bH2 = plot(unqBase2, 5+unqBase2*0);
      set(bH2, 'Color', yColor, 'MarkerFaceColor', yColor);
    end
    bH = [bH1 bH2];
  else
    unqBase = unique(tLaserBaselinePowerMwV);
    if any(unqBase>0);
      bH = plot(unqBase, 5+unqBase*0);
      set(bH, 'Color', 'b', 'MarkerFaceColor', 'b');
    end
  end
  if ~isempty(bH)
    set(bH, ...
        'LineStyle', 'none', ...
        'Marker', 'v', ...
        'MarkerSize', 9);
  end

  anystack([pH1;pH2], 'bottom');
  if all(~isempty(pH))
    anystack(pH, 'top');
  end

  if showLaserStim
    xlabel('power (mW)')
  elseif input.doContrastDetect 
    xlabel('contrast change (%)');
    title(sprintf('base contrast %g', input.baseGratingContrast));
  elseif input.doOriDetect
    xlabel('direction (deg)');
    title(sprintf('base direction %g', input.baseGratingDirectionDeg));
  elseif input.doSpeedDetect
    xlabel('speed (dps)');
    title(sprintf('base speed %g', input.baseGratingSpeedDPS));
  end

  % manually compute limits and tick positions
  xLimL = [min(powerLevels(~powerLevels==0))*0.5, max(powerLevels)*1.5];
  if length(xLimL)==1, xLimL = get(gca, 'XLim'); end
  xL1 = [floor(log10(xLimL(1))) ceil(log10(xLimL(2)))];
  xTickL = 10.^(xL1(1):1:xL1(2));
  xTickL = xTickL(xTickL>=xLimL(1) & xTickL<=xLimL(2));
  if length(xTickL) == 0, 
    % no log10 ticks in range, create two and set them
    xTickL = [xLimL(1) xLimL(2)]; %[xTL(1)/10 xTL(1) xTL(1)]
    xTickL = chop(xTickL,2);
  end
  xTLabelL = cellstr(num2str(xTickL(:)));
  
  set(gca, 'YLim', [0 100], ...
           'XLim', xLimL, ...
           'XScale', 'log', ...
           'XGrid', 'on');
  
  set(gca, 'XTick', xTickL);
  %xTL = chop(get(gca, 'XTick'),2);
  
  set(gca, 'XTickLabel', xTLabelL);
end

%%%%%%%%%%%%%%%%

%% 10 - for several stim levels, plot react times
if nStims > 1
  axH = subplot(spSz{:}, 10);
  hold on;

  powerP.rtMean = [];
  powerP.rtStd = [];

  for iP=1:nL
    tP = powerLevels(iP);
    trIx = stimPowerV == tP;
    
    powerP.rtMean(iP) = mean(reactV(trIx & successIx));
    powerP.rtStd(iP) = std(reactV(trIx & successIx));
    powerP.rtMeanRecent(iP) = mean(reactV(trIx & successIx & lastFewTrIx));
    powerP.rtStdRecent(iP) = std(reactV(trIx & successIx & lastFewTrIx));
    powerP.rtTotalRecent(iP) = sum(trIx & successIx & lastFewTrIx);
    
    if showBlock2
      powerP.rtMean1Block2(iP) = mean(reactV(trIx & successIx & block2Tr1Ix));
      powerP.rtMean2Block2(iP) = mean(reactV(trIx & successIx & block2Tr2Ix));    
      powerP.rtStd1Block2(iP) = std(reactV(trIx & successIx & block2Tr1Ix));
      powerP.rtStd2Block2(iP) = std(reactV(trIx & successIx & block2Tr2Ix));    
      powerP.rtTotal1Block2(iP) = sum(trIx & successIx & block2Tr1Ix);
      powerP.rtTotal2Block2(iP) = sum(trIx & successIx & block2Tr2Ix);    
    end
  end

  % errorbars, manual
  if showBlock2
    des1Ix = ~isnan(powerP.rtMean1Block2);
    des2Ix = ~isnan(powerP.rtMean2Block2);
    h1 = []; h2 = []; hb1 = []; hb2 = [];
    if sum(des1Ix) > 0
      h1 = plot(powerLevels(des1Ix), powerP.rtMean1Block2(des1Ix));
      hb1 = subManualErrorbar(h1, powerP.rtStd1Block2(des1Ix) ./ sqrt(powerP.rtTotal1Block2(des1Ix)));
      set([h1 hb1], 'Color', 'b');
    end
    if sum(des2Ix) > 0 
      h2 = plot(powerLevels(des2Ix), powerP.rtMean2Block2(des2Ix));
      hb2 = subManualErrorbar(h2, powerP.rtStd2Block2(des2Ix) ./ sqrt(powerP.rtTotal2Block2(des2Ix)));
      set([h2 hb2], 'Color', yColor);
    end
    lH = cat(1, h1, h2);
    ebH = cat(1, hb1, hb2);
  else
    lH = plot(powerLevels, powerP.rtMean);
    ebH = subManualErrorbar(lH, powerP.rtStd ./ sqrt(powerP.correct));
  end

  if ~isempty(lH)
    set(lH, 'Marker', 'none', 'LineWidth', 2);
  end

  % recent
  hold on;
  lRH = plot(powerLevels, powerP.rtMeanRecent, '-k');
  ebRH = subManualErrorbar(lRH, powerP.rtStdRecent ./ sqrt(powerP.rtTotalRecent));
  if showBlock2
    set(lRH, 'LineStyle', 'none', ...
             'Marker', 'x');
  end

  % stacking and fixups
  if ~isempty(lRH), 
    anystack(lRH, 'top'); 
  end
  if ~isempty(lH) 
    anystack(lH, 'top'); 
  end
  
  if showLaserStim
    xlabel('power (mW)');
  elseif input.doContrastDetect 
    xlabel('contrast change (%)');
  elseif input.doOriDetect
    xlabel('direction (deg)')
  elseif input.doSpeedDetect
    xlabel('speed (dps)')  
  end
  
  ylabel('rt (ms) on corrects');

  set(gca, ...
      'XLim', xLimL, ...
      'XScale', 'log', ...
      'XGrid', 'on');
  set(gca, 'XTick', xTickL, ...
           'XTickLabel', xTLabelL);  % xTL from above
end

%%%%%%%%%%%%%%%%

%% 12: error rate with hold time?
axH = subplot(spSz{:}, 12);
maxReqX = double(max(input.fixedReqHoldTimeMs)+ ...
          max(input.randReqHoldMaxMs));
if nVisPts > 50
  binWidth = 100; % Mark's code crashed -- pulger
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
clear('pH', 'axH')

[axH] = plotyy(0,1,0,1);
set(axH, 'NextPlot', 'add');

% axes 1
pH(1) = plot(axH(1), edges, N1./totN*100, 'k');
pH(2) = plot(axH(1), edges, N2./totN*100, 'c--');
pH(3) = plot(axH(1), edges, N3./totN*100, 'm--');
set(pH, 'LineWidth', 2)

set(axH(1), 'YLim', [0 100]);
ylabel(axH(1), '% of trials');
title(axH(1), 'trial outcome by hold time');

%axes(axH(2)); hold on;
% axes 2
nPts = 50;
edges=linspace(min(reqHoldV),max(reqHoldV),nPts);
N = histc(reqHoldV, edges);
g = smooth(N, 20, 'lowess');
y0 =  g./max(N);
plot(axH(2), edges, y0, 'b');
set(axH(2), 'YLim', [0 1])

set(axH, 'XLim', [0 ceil(maxReqX/500)*500], ...
         'YTickMode', 'auto');

set(axH(1), 'YColor', 'k');
set(axH(2), 'YColor', 'b');
anystack(axH(2), 'top');



%%%%%%%%%%%%%%%%
%% testing plots


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

%%%%%%%%%%%%%%%%
% debug numbers of trials/odds

%axH = subplot(4,3,9);
%hold on;
%for i=1:8
%    for j=1:2
%        nT(i,j) = sum(tStimulusNumberV == (i-1) & tBlock2TrialNumberV == (j-1));
%    end
%end
%
%plot( (nT./nTrial*80) , '.-');
%set(gca, 'Visible', 'on', ...
%         'YLim', [0 15])
%grid on;

%%%%%%%%%%%%%%%%


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
else
  tV.gratingWidthDeg = input.gratingWidthDeg;
  tV.gratingHeightDeg = input.gratingHeightDeg;
  tV.gratingAzimuthDeg = input.gratingAzimuthDeg;
  tV.gratingElevationDeg = input.gratingElevationDeg;
  tV.gratingSpatialFreqCPD = input.gratingSpatialFreqCPD;
  tV.gratingDurationMs = input.gratingDurationMs;
end

outStr = sprintf('%gx%gdeg, at (%g,%g), %gcpd, %dms', ...
                 roundn(double(tV.gratingWidthDeg),-2), ...
                 roundn(double(tV.gratingHeightDeg),-2), ...
                 tV.gratingAzimuthDeg, ...
                 tV.gratingElevationDeg, ...
                 tV.gratingSpatialFreqCPD, ...
                 tV.gratingDurationMs);


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



