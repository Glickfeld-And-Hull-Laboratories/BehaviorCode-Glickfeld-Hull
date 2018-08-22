function [p, valHR, invHR, valCI, invCI] = allTrialsAttnTest(tValTargets,tInvTargets,hit,miss,invHit,invMiss)

invInd = tInvTargets > 0;
nInvHit = sum(invInd & invHit);
nInvMiss = sum(invInd & invMiss);
[invHR, invCI] = binofit(nInvHit,nInvHit+nInvMiss);

matches = cell2mat(getMatchedValidTrialIndex(tValTargets,tInvTargets(invInd)));

nValHit = sum(hit(matches));
nValMiss = sum(miss(matches));
[valHR,valCI] = binofit(nValHit,nValHit+nValMiss);

p = binocdf(nInvHit,nInvHit+nInvMiss,valHR);
end
