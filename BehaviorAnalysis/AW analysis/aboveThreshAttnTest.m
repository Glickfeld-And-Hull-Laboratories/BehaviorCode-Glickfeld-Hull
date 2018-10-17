function p = aboveThreshAttnTest(thresh,tValTargets,tInvTargets,...
    hit,miss,invHit,invMiss)
    
    tInvAboveThresh = tInvTargets > 0 & tInvTargets > thresh;
    
    nInvHit = sum(tInvAboveThresh & invHit);
    nInvMiss = sum(tInvAboveThresh & invMiss);
    
    matches = cell2mat(getMatchedValidTrialIndex(tValTargets,...
        tInvTargets(tInvAboveThresh)));
    
    nValHit = sum(hit(matches));
    nValMiss = sum(miss(matches));
    valLowHR = nValHit./(nValHit+nValMiss);
    
    p = binocdf(nInvHit,nInvHit+nInvMiss,valLowHR);
    
end
