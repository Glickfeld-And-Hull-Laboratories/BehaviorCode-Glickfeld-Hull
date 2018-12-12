function [p, valHR, invHR, valCI, invCI] = belowThreshAttnTest(thresh,tValTargets,tInvTargets,...
    hit,miss,invHit,invMiss)
    
    tInvBelowThresh = tInvTargets > 0 & tInvTargets <= thresh;
    
    nInvHit = sum(tInvBelowThresh & invHit);
    nInvMiss = sum(tInvBelowThresh & invMiss);
    [invHR, invCI] = binofit(nInvHit,nInvHit+nInvMiss);
    
    matches = cell2mat(getMatchedValidTrialIndex(tValTargets,...
        tInvTargets(tInvBelowThresh)));
    
    nValHit = sum(hit(matches));
    nValMiss = sum(miss(matches));
    valLowHR = nValHit./(nValHit+nValMiss);
    [valHR,valCI] = binofit(nValHit,nValHit+nValMiss);
    
    p = binocdf(nInvHit,nInvHit+nInvMiss,valLowHR);
    
end
