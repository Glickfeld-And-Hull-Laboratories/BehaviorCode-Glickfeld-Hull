function outS = fnComputeThresholdRtMean(varargin)
% computes the mean reaction time at threshold via linear regression through
% RTs at different intensities.
%
%  This is an example use of the data analysis functions
%
% histed 121114


userDefs = { ...
    'Subject', [], ...
    'DateStr', '', ...
    'UniqueId', '', ...
    'LoadDataDebug', false };
    
uo = stropt2struct(stropt_defaults(userDefs, varargin));

bc = behavConstsHADC8;

%% load data
[bs ds x1d fitX1d xd fitXd fitS bootS] ...
    = getAllBehavDataUsingIndex(...
    'Subject', uo.Subject, ...
    'DateStr', uo.DateStr, ...
    'UniqueId', uo.UniqueId, ...
    'LoadDataDebug', uo.LoadDataDebug);


for iBlock = 1:bs.nBlock2Indices
    % exclude low pct corr
    rtPctCorr = feval(fitS(iBlock).modelFun, fitS(iBlock).coefEsts, bs.intensitiesC{iBlock});
    desRtIx = rtPctCorr >= 0.40;
    % linear regression in log space
    logX = log10(bs.intensitiesC{iBlock});
    X = [logX(desRtIx)', ones([sum(desRtIx), 1])];
    Y = bs.reactTimeMean(desRtIx,iBlock);
    [b, bint] = regress(Y, X);
    
    % 68% RT
    tThreshMw = fitX1d.(sprintf('Threshold%d', iBlock));
    rtAtThresh(iBlock) = b(1)*log10(tThreshMw) + b(2);
    
    % find first power above 68% corr, find mean/var
    firstPctCorrN = find(rtPctCorr > 0.68, 1, 'first');
    rtV = bs.reactTimesByPower{firstPctCorrN, iBlock};
    %length(rtV)
    %assert(length(rtV) >= 8, 'short day or bug?');
    
    outS.rtATVar(iBlock) = var(rtV);
    outS.rtATMean(iBlock) = mean(rtV);
end

outS.rtAtThresh = rtAtThresh;

return

%% test this
fnComputeThresholdRtMean('Subject', 48, 'DateStr', '120901', 'UniqueId', '')
