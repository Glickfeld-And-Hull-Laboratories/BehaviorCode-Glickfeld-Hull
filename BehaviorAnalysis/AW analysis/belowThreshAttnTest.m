function p = belowThreshAttnTest(thresh,tValTargets,tInvTargets,...
    hit,miss,invHit,invMiss)
    
    tInvBelowThresh = tInvTargets > 0 & tInvTargets <= thresh;
    
    nInvHit = sum(tInvBelowThresh & invHit);
    nInvMiss = sum(tInvBelowThresh & invMiss);
    
    matches = cell2mat(getMatchedValidTrialIndex(tValTargets,...
        tInvTargets(tInvBelowThresh)));
    
    nValHit = sum(hit(matches));
    nValMiss = sum(miss(matches));
    valLowHR = nValHit./(nValHit+nValMiss);
    
    p = binocdf(nInvHit,nInvHit+nInvMiss,valLowHR);
    
end
