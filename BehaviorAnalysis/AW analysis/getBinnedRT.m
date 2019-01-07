function [RT, RTSte, targetsBinned, targetsBinnedErr, RTanovaP] = ...
    getBinnedRT(tTargets,binEdges,nBins,minTrials,tRT,hit)
    
%     targets = unique(tTargets);
    [~,~,binInd] = histcounts(tTargets,binEdges);
    
%     nHits = nan(1,nBins);
%     nMisses = nan(1,nBins);
    RT = nan(1,nBins);
    RTSte = nan(1,nBins);
    targetsBinned = nan(1,nBins);
    targetsBinnedErr = nan(1,nBins);
    tRT_minTrCorrect = tRT;
    tRT_minTrCorrect(binInd == 0 | ~hit) = nan;
    tRT_allTrialsInBins = cell(1,nBins);
    for ibin = 1:nBins
        ind = hit & binInd == ibin & tRT < 550;
        if sum(ind) > minTrials
            RT(ibin) = mean(tRT(ind));
            RTSte(ibin) = ste(tRT(ind),2);
            targetsBinned(ibin) = mean(tTargets(ind));
            targetsBinnedErr(ibin) = ste(tTargets(ind),2);
            tRT_allTrialsInBins{ibin} = tRT(ind);
        else
            tRT_minTrCorrect(ind) = nan;
            tRT_allTrialsInBins{ibin} = [];
        end
    end
    if nargout > 4
        ind = hit & tTargets > 0;
        RTanovaP = anova1(tRT(ind),binInd(ind));
    end


end