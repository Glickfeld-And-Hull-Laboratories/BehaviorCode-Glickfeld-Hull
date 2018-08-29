function [valRTs_matched, invRTs_matched] = getMatchedReactTimes(...
    tValTargets,tInvTargets,valRT,invRT,valHits,invHits)

valInd = tValTargets > 0 & valHits;
invInd = tInvTargets > 0 & invHits;
matches = cell2mat(getMatchedValidTrialIndex(...
    tValTargets(valInd),tInvTargets(invInd)));

valRT_ind = valRT(valInd);
invRT_ind = invRT(invInd);

valRTs_matched = valRT_ind(matches);
invRTs_matched = invRT_ind;

end