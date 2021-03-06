function [outDs bootDs axH] = weibullFitWithBootstrap(varargin)
%
%   [outDs bootStats] = weibullFitWithBootstrap(varargin)
%
%   Note - this uses *simultaneous* confidence intervals - i.e. for a 95% CI,
%   the range is from 2.5% - 97.5%.
%
% histed 120528 - convert to stropts

%% opts
userDefs = { ...
    'Filename', [], ...
    'BehavDataStruct', [], ... % specify one of Filename, BDS, or X/YPts/NTrials
    'XPoints', [], ...  % specify these three if not reading data from disk
    'YPoints', [], ...
    'NTrials', [], ...
    'DoPlot', false, ...
    'DoPlotAnnotations', true, ...
    'WhichTrials', [], ...
    'DoCorrectEarlies', false, ...
    'DoClampAtZero', true, ...
    'DoBootstrap', false, ...
    'NBootstrapReps', 1000, ...
    'Debug', false, ...
    'DataIndex', [], ...
    'DoExtraBootDistPlots', false};
    
uo = stropt2struct(stropt_defaults(userDefs, varargin));

%% get data if requested
haveBs = false;
if  ~isempty(uo.Filename)
    % get data from mat file
    bs = getBehavDataForDay(...
        'Filename', uo.Filename, ...
        'WhichTrials', uo.WhichTrials, ...
        'DoCorrectEarlies', uo.DoCorrectEarlies, ...
        'DataIndex', uo.DataIndex);
    haveBs = true;
    
    assert(isempty(uo.BehavDataStruct), 'specify a single data source: filename or bds');
elseif ~isempty(uo.BehavDataStruct)
    bs = uo.BehavDataStruct;
    haveBs = true;
    assert(isempty(uo.Filename), 'specify a single data source: filename or bds');
end

if haveBs
    if bs.nBlock2Indices > 1
        error('Only one type of trials should be included in fit');
    end
    
    intensV = bs.intensitiesC{1};
    pctCorr = bs.percentsCorrect;
    nTrialsPerLevel = bs.nCorrPlusMiss;
    
    dropIx = nTrialsPerLevel == 0;
    if any(dropIx)
        intensV = intensV(~dropIx);
        pctCorr = pctCorr(~dropIx);
        nTrialsPerLevel = nTrialsPerLevel(~dropIx);
    end
end


%% supplant with x, y
if ~isempty(uo.XPoints)
    intensV = uo.XPoints;
end
if ~isempty(uo.YPoints)
    pctCorr = uo.YPoints;
end
if ~isempty(uo.NTrials)
    nTrialsPerLevel = uo.NTrials;
end

    

%% run fit
fitS = weibullFit(intensV, pctCorr, uo.DoClampAtZero, false, ...
    { 'nTrials', nTrialsPerLevel'} ); % weight by # trials

slope = fitS.coefEsts(2);

outDs.slope = fitS.slope;
outDs.coefEsts = fitS.coefEsts;
outDs.modelFun = fitS.modelFun;
outDs.fitWeights = fitS.fitWeights;

%% fig
if uo.DoPlot
    if ~ishold; 
        figH = figure; 
        hold on; 
    end
    isContrast = false;
    if exist('bs')
        if bs.intensityIsChr2
            unitStr = 'mW';
            labelStr = 'power (mW)';
            intensScaleF = 1;
        else
            unitStr = '%';
            labelStr = 'contrast (%)';
            intensScaleF = 100;
            isContrast = true;
        end
    else
        unitStr = '';
        labelStr = '';
        intensScaleF = 1;
    end
    
    
    maxI = max(intensV);
    minI = min(intensV); 
    xgrid = logspace(log10(minI*0.5),log10(maxI*1.5),100);
    line(xgrid*intensScaleF, fitS.modelFun(fitS.coefEsts, xgrid), 'Color','k');
    hold on;
    plot(intensV*intensScaleF, pctCorr, 'o');
    %thresh = coefEsts(1)*[1 1];
    plot(fitS.thresh*[1 1]*intensScaleF, [0 fitS.threshY], '--');
    
    % set limits correctly
    xLim = [min(xgrid) max(xgrid)] .* intensScaleF .* [0.75 1.25];
    xLim = 10.^ceil(log10(xLim) - [1 0]);
    
    if isContrast
        xLim(2) = min(xLim(2), 100);
    end
    set(gca, 'XLim', xLim);
    set(gca, 'XScale', 'log');
    set(gca, 'YLim', [0 1]);
    
    if uo.DoPlotAnnotations
        tStr = sprintf(['Threshold: %3g %s\n', ...
            'Slope: %3g'], ...
            chop(fitS.thresh,2)*intensScaleF, ...
            unitStr, ...
            chop(fitS.slope,2));
        text(xLim(2), 0.15, tStr, ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');
    end
    
    xTick = get(gca, 'XTick');
    xTL = eval(['{' sprintf('''%g'',', chop(xTick,2)) '}']);
    set(gca, 'XTickLabel', xTL);
    xlabel(labelStr);
    ylabel('fraction correct');

end

outDs = struct_union(outDs, fitS);
if ~uo.DoBootstrap
    bootStats = []; 
end

%return

%% do bootstrap
if uo.DoBootstrap
    if uo.DoPlot
        drawnow;
    end
    if uo.Debug
        disp(sprintf('Running %d bootstrap reps', uo.NBootstrapReps));
    end
    
    % first pull out nCorr and nTot for all
    nTotB = nTrialsPerLevel;

    pctCorrB = pctCorr;
    pctCorrB(pctCorrB>=1) = 1-10*eps;
    pctCorrB(pctCorrB<=0) = 0+10*eps;
    
    threshB = repmat(NaN, [uo.NBootstrapReps, 1]);
    coefEstsB = repmat(NaN, [uo.NBootstrapReps, 4]);
    weightsB = repmat(NaN, [uo.NBootstrapReps length(nTotB)]);
    
    for iR = 1:uo.NBootstrapReps
        tNCorr = binornd(nTotB, pctCorrB);
    
        tPctCorr = tNCorr./nTotB;
        tPctCorr(tPctCorr >= 1) = 1-10*eps;
        tPctCorr(tPctCorr == 0) = 0+10*eps;
        
        % do the fit
        bootFitS = weibullFit(intensV, tPctCorr, uo.DoClampAtZero, false, ...
            {'nTrials', nTotB });  
        % 120710 - I choose to allow the weights to be recalculated from the
        % pctCorr on each bootstrap trial.  It probably yields tighter
        % intervals to fix the weights based on the actual data sample.

        threshB(iR) = bootFitS.thresh;
        coefEstsB(iR,:) = bootFitS.coefEsts;
        weightsB(iR,:) = bootFitS.fitWeights(:)';
        
        if mod(iR, 25) == 0
            statusbar(0, '%s done %4g%%: %d of %d bootstrap reps', ...
                mfilename, iR/uo.NBootstrapReps*100, iR, uo.NBootstrapReps);
        end
    end
    
    nNaN = mean(sum(isnan(coefEstsB)));
    if nNaN > 0
        disp(sprintf('** Errors/Reps: %d/%d ', nNaN, uo.NBootstrapReps));
    end
    
    slopeB = coefEstsB(:,2);

    if nargout > 1
        bootDs.coefEstsMat = coefEstsB;
        bootDs.fitWeightsMat = weightsB;
    end
    %keyboard

    bootStats.threshMean = nanmean(threshB);  % should be near the true value
    bootStats.threshStd = nanstd(threshB);  % should be near the true value
    bootStats.ci95 = prctile(threshB, [0 100] + 5*[1 -1]/2);
    bootStats.ci99 = prctile(threshB, [0 100] + 1*[1 -1]/2);

    bootStats.slopeMean = nanmean(slopeB);
    bootStats.slopeStd = std(slopeB);
    bootStats.slopeCi95 = prctile(slopeB, [0 100] + 5*[1 -1]/2);
    bootStats.slopeCi99 = prctile(slopeB, [0 100] + 1*[1 -1]/2);
    
    if uo.DoPlot
        plot(bootStats.ci95*intensScaleF, fitS.threshY*[1 1], 'k-');
        if nNaN > 0
            text(xLim(2), 0.5, ...
                sprintf('Num bootstrap fit errors: %d/%d', nNaN, uo.NBootstrapReps), ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'right');
        end
        
        if uo.DoPlotAnnotations
            tStr = sprintf('95%% CI, thresh: %3g-%3g\n95%% CI, slope: %3g-%3g', ...
                chop(bootStats.ci95(1), 2)*intensScaleF, ...
                chop(bootStats.ci95(2), 2)*intensScaleF, ...
                chop(bootStats.slopeCi95(1),2), ...
                chop(bootStats.slopeCi95(2),2));
            text(xLim(2), 0.05, tStr, ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'right');
        end
    end

    outDs.bootStats = bootStats;
        
    if uo.DoExtraBootDistPlots
        figH = figure;
        subplot(1,2,1);
        cdfplot(threshB);
        xlabel('threshold')
        subplot(1,2,2);
        cdfplot(slopeB);
        xlabel('slope');
        keyboard
    end
end





