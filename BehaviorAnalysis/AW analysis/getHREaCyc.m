function [valHR,invHR] = getHREaCyc(msCmlvData,valTargets,invTargets,maxCycles,minTrN_timebin)

    ncyc = maxCycles;

    cycles = unique(cat(2,msCmlvData.valTargetCycle(~isnan(msCmlvData.valTargetCycle)),...
        msCmlvData.invTargetCycle(~isnan(msCmlvData.invTargetCycle)))+1);
    
%     valTimeID = discretize(msCmlvData.valTargetTimeMs,targetTimeBins);
%     invTimeID = discretize(msCmlvData.invTargetTimeMs,targetTimeBins);

    valHR = nan(1,ncyc);
    invHR = nan(1,ncyc);
    for icyc = cycles
        if icyc > maxCycles
            continue
        elseif icyc == maxCycles
            ind1 = msCmlvData.valTargetCycle+1 >= icyc;    
            ind2 = msCmlvData.invTargetCycle+1 >= icyc;    
        else
            ind1 = msCmlvData.valTargetCycle+1 == icyc;    
            ind2 = msCmlvData.invTargetCycle+1 == icyc;     
        end
        invInd = ind2 & invTargets > 0 & (msCmlvData.invHit|msCmlvData.invMiss);
        matches = cell2mat(getMatchedValidTrialIndex(valTargets,...
            invTargets(invInd)));
        valInd = matches(ismember(matches,find(ind1)));
        if sum(ind1) > minTrN_timebin & sum(ind2) > minTrN_timebin 
            valHR(icyc) = sum(msCmlvData.hit(valInd))./...
                sum(msCmlvData.hit(valInd)|msCmlvData.miss(valInd));
            invHR(icyc) = sum(msCmlvData.invHit(invInd))./...
                sum(msCmlvData.invHit(invInd)|msCmlvData.invMiss(invInd));
        end
    end
end