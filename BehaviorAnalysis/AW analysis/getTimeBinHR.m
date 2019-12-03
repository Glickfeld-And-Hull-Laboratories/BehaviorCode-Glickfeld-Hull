function [valHR,invHR] = getTimeBinHR(msCmlvData,valTargets,invTargets,targetTimeBins)

    nTimeBins = length(targetTimeBins)-1;

    valTimeID = discretize(msCmlvData.valTargetTimeMs,targetTimeBins);
    invTimeID = discretize(msCmlvData.invTargetTimeMs,targetTimeBins);

    valHR = nan(1,nTimeBins);
    invHR = nan(1,nTimeBins);
    for ibin = 1:nTimeBins
        ind1 = valTimeID == ibin;    
        ind2 = invTimeID == ibin;     
        invInd = ind2 & invTargets > 0 & (msCmlvData.invHit|msCmlvData.invMiss); 
        matches = cell2mat(getMatchedValidTrialIndex(valTargets,...
            invTargets(invInd)));
        valInd = matches(ismember(matches,find(ind1)));
        valHR(ibin) = sum(msCmlvData.hit(valInd))./...
            sum(msCmlvData.hit(valInd)|msCmlvData.miss(valInd));
        invHR(ibin) = sum(msCmlvData.invHit(invInd))./...
            sum(msCmlvData.invHit(invInd)|msCmlvData.invMiss(invInd));
    end
end