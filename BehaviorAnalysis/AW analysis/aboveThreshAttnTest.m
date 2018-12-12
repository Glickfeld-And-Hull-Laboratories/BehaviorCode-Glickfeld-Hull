function [p, valHR, invHR, valCI, invCI] = aboveThreshAttnTest(thresh,tValTargets,tInvTargets,...
    hit,miss,invHit,invMiss)
    
    tInvAboveThresh = tInvTargets > 0 & tInvTargets > thresh;
    
    nInvHit = sum(tInvAboveThresh & invHit);
    nInvMiss = sum(tInvAboveThresh & invMiss);
    [invHR, invCI] = binofit(nInvHit,nInvHit+nInvMiss);
    
    matches = cell2mat(getMatchedValidTrialIndex(tValTargets,...
        tInvTargets(tInvAboveThresh)));
    
    nValHit = sum(hit(matches));
    nValMiss = sum(miss(matches));
    valLowHR = nValHit./(nValHit+nValMiss);
    [valHR,valCI] = binofit(nValHit,nValHit+nValMiss);
    
    p = binocdf(nInvHit,nInvHit+nInvMiss,valLowHR);
    
end
