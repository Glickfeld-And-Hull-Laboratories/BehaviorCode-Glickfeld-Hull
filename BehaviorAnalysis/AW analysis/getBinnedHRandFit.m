function [HR, HRci, targetsBinned, targetsBinnedErr, HRfit] = ...
    getBinnedHRandFit(tTargets,binEdges,nBins,minTrials,hit,miss,FAR,nFAR)
    
%     targets = unique(tTargets);
    [~,~,binInd] = histcounts(tTargets,binEdges);
    
    nHits = nan(1,nBins);
    nMisses = nan(1,nBins);
    targetsBinned = nan(1,nBins);
    targetsBinnedErr = nan(1,nBins);
    for ibin = 1:nBins
        ind = (hit | miss) & binInd == ibin;
        if sum(ind) > minTrials
            nHits(ibin) = sum(hit & binInd == ibin);
            nMisses(ibin) = sum(miss & binInd == ibin);
            targetsBinned(ibin) = mean(tTargets(ind));
            targetsBinnedErr(ibin) = ste(tTargets(ind),2);
        end
    end
    nanInd = ~isnan(nHits);
    targetsBinned = targetsBinned(nanInd);
    targetsBinnedErr = targetsBinnedErr(nanInd);
    
    if sum(nHits(nanInd)+nMisses(nanInd)) > 0
        [HR,HRci] = binofit(nHits(nanInd),nHits(nanInd)+nMisses(nanInd));
    else
        HR = nan;
        HRci = nan(1,2);
    end
   
    if nargout == 5
        maxI = max(targetsBinned);
        minI = min(targetsBinned);
    %     xGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
        if isempty(FAR)
            nTrials = nHits(nanInd)+nMisses(nanInd);
            HRfit = weibullFitLG(targetsBinned, HR, 0,0, {'nTrials',nTrials});
        else
            nTrials = [nFAR,nHits(nanInd)+nMisses(nanInd)];
            HRfit = weibullFitLG([0, targetsBinned], [FAR, HR], 0 ,0, {'nTrials',nTrials});
        end
    end
end