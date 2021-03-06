function [valHR_highThreshold, invHR_highThreshold,...
    valCI_highThreshold, invCI_highThreshold, nVal, nInv, ...
    msHighThresh,val_matchTrialHiLoIndex] = ...
    getMatchedHighThresholdHR(targets,HRFit,highThreshold,valTargets,invTargets,...
    hit,miss,invHit,invMiss)

    maxI = max(targets);
    minI = min(targets);
    fitGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    highThreshInd = find(HRFit.modelFun(HRFit.coefEsts,fitGrid) ...
        > highThreshold,1);
    if isempty(highThreshInd)
        highThreshInd = find(HRFit.modelFun(HRFit.coefEsts,fitGrid) ...
            > 0.85,1);
    end
    msHighThresh = fitGrid(highThreshInd);

    invHR_highThreshold= nan(1,2);
    t = invTargets < msHighThresh & invTargets > 0;
    invHR_highThreshold(1) = sum(t & invHit)./sum(t & (invHit | invMiss));
    nInv(1) = sum(t & (invHit | invMiss));
    [~,invCI_highThreshold(1,:)] = binofit(sum(t & invHit),sum(t & (invHit | invMiss)));
    t = invTargets >= msHighThresh;
    invHR_highThreshold(2) = sum(t & invHit)./sum(t & (invHit | invMiss));
    nInv(2) = sum(t & (invHit | invMiss));
    [~,invCI_highThreshold(2,:)] = binofit(sum(t & invHit),sum(t & (invHit | invMiss)));
    invCI_highThreshold = invCI_highThreshold';
    
    matchedTargets = ismember(valTargets,unique(invTargets));

    valHR_highThreshold = nan(1,2);
    valCI_highThreshold = nan(2,2);
    t = valTargets < msHighThresh & valTargets > 0 & matchedTargets;
    val_matchTrialHiLoIndex = zeros(1,length(valTargets));
    val_matchTrialHiLoIndex(t) = 1;
    valHR_highThreshold(1) = sum(t & hit)./sum(t & (hit | miss));
    nVal(1) = sum(t & (hit | miss));
    [~,valCI_highThreshold(1,:)] = binofit(sum(t & hit),sum(t & (hit | miss)));
    t = valTargets >= msHighThresh & matchedTargets;
    valHR_highThreshold(2) = sum(t & hit)./sum(t & (hit | miss));
    nVal(2) = sum(t & (hit | miss));
    [~,valCI_highThreshold(2,:)] = binofit(sum(t & hit),sum(t & (hit | miss)));
    val_matchTrialHiLoIndex(t) = 2;
    valCI_highThreshold = valCI_highThreshold';
    
        
end
