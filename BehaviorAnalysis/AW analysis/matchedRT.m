function [valRT,invHR,valSEM,invSEM] = matchedRT(tValTargets,tInvTargets,...
    tValRT,tInvRT,invHit,invMiss)

    invInd = tInvTargets > 0;
    invRT = mean(tInvRT(invInd));
    invSEM = ste(tInvRT(invInd),1);
    
    matches = cell2mat(getMatchedValidTrialIndex(tValTargets,tInvTargets(invInd)));
    
    valRT = mean(tValRT(matches));
    valSEM = ste(tInvRT(matches),1);
end