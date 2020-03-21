function [valHR,invHR] = getHREaCyc_av(msCmlvData,maxCycles,minTrN_timebin,visBinEdges,audBinEdges)
    rng(0)
    ncyc = maxCycles;
    nBoot = 100;

    cycles = unique(cat(2,msCmlvData.valTargetCycle(~isnan(msCmlvData.valTargetCycle)),...
        msCmlvData.invTargetCycle(~isnan(msCmlvData.invTargetCycle)))+1);
    
%     valTimeID = discretize(msCmlvData.valTargetTimeMs,targetTimeBins);
%     invTimeID = discretize(msCmlvData.invTargetTimeMs,targetTimeBins);

    valTargets_vis = discretize(msCmlvData.tVisTargets,[0 visBinEdges])-1;
    invTargets_vis = discretize(msCmlvData.tInvVisTargets,[0 visBinEdges])-1;
    valTargets_aud = discretize(msCmlvData.tAudTargets,[0 audBinEdges])-1;
    invTargets_aud = discretize(msCmlvData.tInvAudTargets,[0 audBinEdges])-1;
    visTargets = unique(invTargets_vis(~isnan(invTargets_vis)&invTargets_vis~=0));
    audTargets = unique(invTargets_aud(~isnan(invTargets_aud)&invTargets_aud~=0));
    
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
        if sum(ind1)<=minTrN_timebin | sum(ind2)<=minTrN_timebin
            continue
        end
        
        invInd_vis = ind2 & invTargets_vis > 0 & (msCmlvData.invHit|msCmlvData.invMiss);
        invInd_aud = ind2 & invTargets_aud > 0 & (msCmlvData.invHit|msCmlvData.invMiss);
        
        val_matches = [];
        inv_matches = [];
        for i = visTargets
            if isempty(find(invInd_vis & invTargets_vis == i))
                continue
            elseif sum(invInd_vis & invTargets_vis == i) == 1
                inv_matches = cat(2,inv_matches,repmat(find(invInd_vis & invTargets_vis == i),[1,nBoot]));
                val_matches = cat(2,val_matches,randsample(find(ind1 & valTargets_vis == i & (msCmlvData.hit|msCmlvData.miss)),nBoot,1));
            else
                inv_matches = cat(2,inv_matches,randsample(find(invInd_vis & invTargets_vis == i),nBoot,1));
                val_matches = cat(2,val_matches,randsample(find(ind1 & valTargets_vis == i & (msCmlvData.hit|msCmlvData.miss)),nBoot,1));
            end
        end
        for i = audTargets
            if isempty(find(invInd_aud & invTargets_aud == i))
                continue
            elseif sum(invInd_aud & invTargets_aud == i) == 1
                inv_matches = cat(2,inv_matches,repmat(find(invInd_aud & invTargets_aud == i),[1,nBoot]));
                val_matches = cat(2,val_matches,randsample(find(ind1 & valTargets_aud == i & (msCmlvData.hit|msCmlvData.miss)),nBoot,1));
            else
                inv_matches = cat(2,inv_matches,randsample(find(invInd_aud & invTargets_aud == i),nBoot,1));
                val_matches = cat(2,val_matches,randsample(find(ind1 & valTargets_aud == i & (msCmlvData.hit|msCmlvData.miss)),nBoot,1));
            end
        end
       
%         if sum(ind1) > minTrN_timebin && sum(ind2) > minTrN_timebin 
        valHR(icyc) = sum(msCmlvData.hit(val_matches))./...
            sum(msCmlvData.hit(val_matches)|msCmlvData.miss(val_matches));
        invHR(icyc) = sum(msCmlvData.invHit(inv_matches))./...
            sum(msCmlvData.invHit(inv_matches)|msCmlvData.invMiss(inv_matches));
%         end
    end
end