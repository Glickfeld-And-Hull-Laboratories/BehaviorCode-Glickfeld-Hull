function msSum = getEaMouseBxSummary(ms,nBins,minTrN_ms)
    visualTrials = 1;
    auditoryTrials = 2;
    val = 1;
    inv = 2;
    oriBin = 1;
    timeBin = 2;
    msSummary = struct;
    nmice = size(ms,2);
    longTrialCutoff = 5;
for im = 1:nmice
    tVisTargets = round(double(ms(im).visTargets), 2, 'significant');
    visTrials = tVisTargets > 0;
    tAudTargets = round(ms(im).audTargets, 2, 'significant');
    audTrials = tAudTargets > 0;
    visTargets = unique(tVisTargets);
    audTargets = unique(tAudTargets);
    visTargets = visTargets(2:end);
    audTargets = audTargets(2:end);
    visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
    audBinEdges = exp(linspace(log(min(audTargets)-(0.5*min(audTargets))),log(max(audTargets)),nBins+1));
    [~,~,visBinInd] = histcounts(tVisTargets,visBinEdges);
    [~,~,audBinInd] = histcounts(tAudTargets,audBinEdges);
    hit = ms(im).successIx;
    miss = ms(im).missedIx;

    visHitRate = nan(1,nBins);
    audHitRate = nan(1,nBins);
    nVisHits = nan(1,nBins);
    nVisMisses = nan(1,nBins);
    nAudHits = nan(1,nBins);
    nAudMisses = nan(1,nBins);
    visTargetsBinned = nan(1,nBins);
    audTargetsBinned = nan(1,nBins);
    visTargetsSte = nan(1,nBins);
    audTargetsSte = nan(1,nBins);
    for i = 1:nBins
        ind = (hit | miss) & visBinInd == i;
        nVisHits(i) = sum(hit & visBinInd == i);
        nVisMisses(i) = sum(miss & visBinInd == i);
        if sum(ind) < minTrN_ms
            visTargetsBinned(i) = nan;
            visTargetsSte(i) = nan;
            visHitRate(i) = nan;
        else
            visTargetsBinned(i) = mean(tVisTargets(ind));
            visTargetsSte(i) = ste(tVisTargets(ind),2);
            visHitRate(i) = sum(hit(visBinInd == i))./ ...
                sum(hit(visBinInd == i) | miss(visBinInd == i));
        end
        ind = (hit | miss) & audBinInd == i;
        nAudHits(i) = sum(hit & audBinInd == i);
        nAudMisses(i) = sum(miss & audBinInd == i);
        if sum(ind) < minTrN_ms
            audTargetsBinned(i) = nan;
            audTargetsSte(i) = nan;
            audHitRate(i) = nan;
        else
            audTargetsBinned(i) = mean(tAudTargets(ind));
            audTargetsSte(i) = ste(tAudTargets(ind),2);
            audHitRate(i) = sum(hit(audBinInd == i))...
                ./ sum(hit(audBinInd == i) | miss(audBinInd == i));
        end
    end
    [~,visHR95ci] = binofit(nVisHits,nVisHits+nVisMisses);
    [~,audHR95ci] = binofit(nAudHits,nAudHits+nAudMisses);

    tInvVisTargets = round(ms(im).invVisTargets,2,'significant');
    tInvVisTargets(isnan(tInvVisTargets)) = 0;
    invVisTrials = tInvVisTargets > 0;
    tInvAudTargets = round(ms(im).invAudTargets,2,'significant');
    tInvAudTargets(isnan(tInvAudTargets)) = 0;
    invAudTrials = tInvAudTargets > 0;
    [~,~,invVisBinInd] = histcounts(tInvVisTargets,visBinEdges);
    [~,~,invAudBinInd] = histcounts(tInvAudTargets,audBinEdges);
    invHit = ms(im).invHitIx;
    invMiss = ms(im).invMissIx;

    invVisHitRate = nan(1,nBins);
    invAudHitRate = nan(1,nBins);
    nInvVisHits = nan(1,nBins);
    nInvVisMisses = nan(1,nBins);
    nInvAudHits = nan(1,nBins);
    nInvAudMisses = nan(1,nBins);
    invVisTargetsBinned = nan(1,nBins);
    invAudTargetsBinned = nan(1,nBins);
    invVisTargetsSte = nan(1,nBins);
    invAudTargetsSte = nan(1,nBins);
    for i = 1:nBins
        ind = (invHit | invMiss) & invVisBinInd == i;
        nInvVisHits(i) = sum(invVisBinInd == i & invHit);
        nInvVisMisses(i) = sum(invVisBinInd == i & invMiss);
        if sum(ind) < minTrN_ms
            invVisTargetsBinned(i) = nan;
            invVisTargetsSte(i) = nan;
            invVisHitRate(i) = nan;
        else
            invVisTargetsBinned(i) = mean(tInvVisTargets(ind));
            invVisTargetsSte(i) = ste(tInvVisTargets(ind),2);
            invVisHitRate(i) = sum(invHit(invVisBinInd == i))...
                ./sum(invHit(invVisBinInd == i) | invMiss(invVisBinInd == i));
        end
        ind = (invHit | invMiss) & invAudBinInd == i;
        nInvAudHits(i) = sum(invAudBinInd == i & invHit);
        nInvAudMisses(i) = sum(invAudBinInd == i & invMiss);
        if sum(ind) < minTrN_ms
            invAudTargetsBinned(i) = nan;
            invAudTargetsSte(i) = nan;
            invAudHitRate(i) = nan;
        else
            invAudTargetsBinned(i) = mean(tInvAudTargets(ind));
            invAudTargetsSte(i) = ste(tInvAudTargets(ind),2);
            invAudHitRate(i) = sum(invHit(invAudBinInd == i))...
                ./sum(invHit(invAudBinInd == i) | invMiss(invAudBinInd == i));
        end
    end
    [~,invVisHR95ci] = binofit(nInvVisHits,nInvVisHits+nInvVisMisses);
    [~,invAudHR95ci] = binofit(nInvAudHits,nInvAudHits+nInvAudMisses);
    
    %%
    ind = ~isnan(visHitRate);
    msFit = weibullFitLG(visTargetsBinned(ind), visHitRate(ind),...
        1, 0, {'nTrials', nVisHits(ind)+nVisMisses(ind)});

    visBTAttnTestP = belowThreshAttnTest(msFit.thresh,tVisTargets,tInvVisTargets,...
        hit,miss,invHit,invMiss);
    visAttnTestP = allTrialsAttnTest(tVisTargets,tInvVisTargets,...
        hit,miss,invHit,invMiss);
    
    ind = ~isnan(audHitRate);
    msFit = weibullFitLG(audTargetsBinned(ind), audHitRate(ind),...
        1, 0, {'nTrials', nAudHits(ind)+nAudMisses(ind)});

    audBTAttnTestP = belowThreshAttnTest(msFit.thresh,tAudTargets,tInvAudTargets,...
        hit,miss,invHit,invMiss);
    audAttnTestP = allTrialsAttnTest(tAudTargets,tInvAudTargets,...
        hit,miss,invHit,invMiss);
    %%
    %%
    visMatches = cell2mat(getMatchedValidTrialIndex(tVisTargets,tInvVisTargets));
    matchedVisHR = sum(ismember(visMatches,find(hit))) ./ ...
        (sum(ismember(visMatches,find(hit)))+sum(ismember(visMatches,find(miss))));
    matchedInvVisHR = sum(invHit & invVisTrials)./sum(invVisTrials & (invHit |invMiss));
    audMatches = cell2mat(getMatchedValidTrialIndex(tAudTargets,tInvAudTargets));
    matchedAudHR = sum(ismember(audMatches,find(hit))) ./ ...
        (sum(ismember(audMatches,find(hit)))+sum(ismember(audMatches,find(miss))));
    matchedInvAudHR = sum(invHit & invAudTrials)./sum(invAudTrials & (invHit |invMiss));
    
    %%
    tCycles = ms(im).nCycles;
    cycles = unique(tCycles);
    nCycles = max(cycles);
    longTrials = tCycles >= longTrialCutoff;
    
    visCycHR = nan(1,nCycles);
    nVisCycHits = nan(1,nCycles);
    nVisCycMisses = nan(1,nCycles);
    audCycHR = nan(1,nCycles);
    nAudCycHits  = nan(1,nCycles);
    nAudCycMisses = nan(1,nCycles);
    
    invCycles = ms(im).invCycles;
    longInvTrials = invCycles >= longTrialCutoff;
    
    invVisCycHR = nan(1,nCycles);
    nInvVisCycHits = nan(1,nCycles);
    nInvVisCycMisses = nan(1,nCycles);
    invAudCycHR = nan(1,nCycles);
    nInvAudCycHits  = nan(1,nCycles);
    nInvAudCycMisses = nan(1,nCycles);
    
    for icyc = 1:nCycles
        nVisCycHits(icyc) = sum(visTrials & tCycles == icyc & hit);
        nVisCycMisses(icyc) = sum(visTrials & tCycles == icyc & miss);
        ind = visTrials & tCycles == icyc & (hit | miss);
        if sum(ind) > minTrN_ms
            visCycHR(icyc) = sum(visTrials & tCycles == icyc & hit) ./ ...
                sum(visTrials & tCycles == icyc & (hit | miss));
        end
        nAudCycHits(icyc) = sum(audTrials & tCycles == icyc & hit);
        nAudCycMisses(icyc) = sum(audTrials & tCycles == icyc & miss);
        ind = audTrials & tCycles == icyc & (hit | miss);
        if sum(ind) > minTrN_ms
            audCycHR(icyc) = sum(audTrials & tCycles == icyc & hit) ./ ...
                sum(audTrials & tCycles == icyc & (hit | miss));
        end
        
        nInvVisCycHits(icyc) = sum(invVisTrials & tCycles == icyc & invHit);
        nInvVisCycMisses(icyc) = sum(invVisTrials & tCycles == icyc & invMiss);
        ind = invVisTrials & tCycles == icyc & (invHit | invMiss);
        if sum(ind) > minTrN_ms
            invVisCycHR(icyc) = sum(invVisTrials & tCycles == icyc & invHit) ./ ...
                sum(invVisTrials & tCycles == icyc & (invHit | invMiss));
        end
        nInvAudCycHits(icyc) = sum(invAudTrials & tCycles == icyc & invHit);
        nInvAudCycMisses(icyc) = sum(invAudTrials & tCycles == icyc & invMiss);
        ind = invAudTrials & tCycles == icyc & (invHit | invMiss);
        if sum(ind) > minTrN_ms
            invAudCycHR(icyc) = sum(invAudTrials & tCycles == icyc & invHit) ./ ...
                sum(invAudTrials & tCycles == icyc & (invHit | invMiss));
        end
    end
    [~,visCyc95ci] = binofit(nVisCycHits,nVisCycHits+nVisCycMisses);
    [~,audCyc95ci] = binofit(nAudCycHits,nAudCycHits+nAudCycMisses);
    [~,invVisCyc95ci] = binofit(nInvVisCycHits,nInvVisCycHits+nInvVisCycMisses);
    [~,invAudCyc95ci] = binofit(nInvAudCycHits,nInvAudCycHits+nInvAudCycMisses);

    %%
    hitLong = hit(longTrials);
    missLong = miss(longTrials);
    invVisLongTrials = tInvVisTargets > 0 & longInvTrials;
    invAudLongTrials = tInvAudTargets > 0 & longInvTrials;
    visMatches = cell2mat(getMatchedValidTrialIndex(tVisTargets(longTrials),...
        tInvVisTargets(longInvTrials)));
    audMatches = cell2mat(getMatchedValidTrialIndex(tAudTargets(longTrials),...
        tInvAudTargets(longInvTrials)));
    
    matchedVisLongHR = sum(ismember(visMatches,find(hitLong)))./...
        (sum(ismember(visMatches,find(hitLong))) + ...
        sum(ismember(visMatches,find(missLong))));
    matchedInvVisLongHR = sum(invHit & invVisLongTrials)./...
        sum(invVisLongTrials & (invHit | invMiss));
    matchedAudLongHR = sum(ismember(audMatches,find(hitLong)))./...
        (sum(ismember(audMatches,find(hitLong))) + ...
        sum(ismember(visMatches,find(missLong))));
    matchedInvAudLongHR = sum(invHit & invAudLongTrials)./...
        sum(invAudLongTrials & (invHit | invMiss));
    
    %%
    msSummary(im).av(visualTrials).belowThreshAttnTestP = visBTAttnTestP;
    msSummary(im).av(auditoryTrials).belowThreshAttnTestP = audBTAttnTestP;
    msSummary(im).av(visualTrials).allTrialsAttnTestP = visAttnTestP;
    msSummary(im).av(auditoryTrials).allTrialsAttnTestP = audAttnTestP;
    
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).hitRate = ...
        visHitRate;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).hitRateCI = ...
        visHR95ci;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).binTarget = ...
        visTargetsBinned;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).binTargetSte = ...
        visTargetsSte;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).nTrials = ...
        nVisHits+nVisMisses;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).matchedHR = ...
        matchedVisHR;
    msSummary(im).av(visualTrials).binning(oriBin).cue(val).matchedLongHR = ...
        matchedVisLongHR;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).hitRate = ...
        invVisHitRate;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).hitRateCI = ...
        invVisHR95ci;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).binTarget = ...
        invVisTargetsBinned;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).binTargetSte = ...
        invVisTargetsSte;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).nTrials = ...
        nInvVisHits+nInvVisMisses;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).matchedHR = ...
        matchedInvVisHR;
    msSummary(im).av(visualTrials).binning(oriBin).cue(inv).matchedLongHR = ...
        matchedInvVisLongHR;

    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).hitRate = ...
        audHitRate;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).hitRateCI = ...
        audHR95ci;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).binTarget = ...
        audTargetsBinned;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).binTargetSte = ...
        audTargetsSte;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).nTrials = ...
        nAudHits+nAudMisses;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).matchedHR = ...
        matchedAudHR;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).matchedLongHR = ...
        matchedAudLongHR;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).hitRate = ...
        invAudHitRate;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).hitRateCI = ...
        invAudHR95ci;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).binTarget = ...
        invAudTargetsBinned;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).binTargetSte = ...
        invAudTargetsSte;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).nTrials = ...
        nInvAudHits+nInvAudMisses;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).matchedHR = ...
        matchedInvAudHR;
    msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).matchedLongHR = ...
        matchedInvAudLongHR;
    
    msSummary(im).av(visualTrials).binning(timeBin).cue(val).hitRate = ...
        visCycHR;
    msSummary(im).av(visualTrials).binning(timeBin).cue(val).hitRateCI = ...
        visCyc95ci;
    msSummary(im).av(visualTrials).binning(timeBin).cue(val).nTrials = ...
        nVisCycHits+nVisCycMisses;
    msSummary(im).av(visualTrials).binning(timeBin).cue(inv).hitRate = ...
        invVisCycHR;
    msSummary(im).av(visualTrials).binning(timeBin).cue(inv).hitRateCI = ...
        invVisCyc95ci;
    msSummary(im).av(visualTrials).binning(timeBin).cue(inv).nTrials = ...
        nInvVisCycHits+nInvVisCycMisses;
    
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(val).hitRate = ...
        audCycHR;
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(val).hitRateCI = ...
        audCyc95ci;
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(val).nTrials = ...
        nAudCycHits+nVisCycMisses;
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(inv).hitRate = ...
        invAudCycHR;
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(inv).hitRateCI = ...
        invAudCyc95ci;
    msSummary(im).av(auditoryTrials).binning(timeBin).cue(inv).nTrials = ...
        nInvAudCycHits+nInvAudCycMisses;
end
    msSum = msSummary;
end