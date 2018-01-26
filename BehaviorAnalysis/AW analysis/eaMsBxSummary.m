clear all
close all

ms2analyze = {'747'};
nmice = length(ms2analyze);

early_cutoff = 0.5;
lapse_cutoff = 0.85;
minTrN_ms = 5;
nBins = 6;
minCyclesFA = 2;

visualTrials = 1;
auditoryTrials = 2;
valid  = 1;
invalid = 2;

rc = behavConstsAV;

%%
for im = 1:nmice
    mouseName = ms2analyze{im};

    bxDataInfo = table2struct(readtable(rc.eaMouseIndexFilename,'Sheet',mouseName));
    % fnin = rc.behavData;
    fnout = fullfile(rc.ashleyAnalysis,mouseName,'behavior');
    fileName = sprintf('%s_bxSummary_%s',mouseName,datestr(now,'yymmdd'));
    if ~exist(fnout,'dir')
        mkdir(fnout)
    end

    msExptInfo = bxFrmXLS(mouseName,bxDataInfo);
    %%
    nexp = size(msExptInfo,2);
    msCmlvData = struct;
    exptThresh = nan(2,nexp);
    valHRexpt = cell(1,nexp);
    invHRexpt = cell(1,nexp);
    FAR = nan(2,nexp);
    
    nExpAnalyzed = {0;0};
    exptInd = false(1,nexp);
    for iexp = 1:nexp
        if length(msExptInfo(iexp).visHR) > 1
            visHRCutoffPass = msExptInfo(iexp).visHR(end-1:end) < lapse_cutoff;
        else
            visHRCutoffPass = msExptInfo(iexp).visHR;
        end
        if length(msExptInfo(iexp).audHR) > 1
            audHRCutoffPass = msExptInfo(iexp).audHR(end-1:end) < lapse_cutoff;
        else
            audHRCutoffPass = msExptInfo(iexp).audHR < lapse_cutoff;
        end
        if isnan(msExptInfo(iexp).invType)
            continue
        elseif msExptInfo(iexp).pctEarly > early_cutoff | ...
                all(visHRCutoffPass) | all(audHRCutoffPass)
            continue
        end
        exptInd(iexp) = true;

        if strcmp(msExptInfo(iexp).invType,'vis')
            nExpAnalyzed{visualTrials} = nExpAnalyzed{visualTrials}+1;
        else
            nExpAnalyzed{auditoryTrials} = nExpAnalyzed{auditoryTrials}+1;
        end

        tVisTargets = round(double(msExptInfo(iexp).tVisTargets),2,'significant');
        visTrials = tVisTargets > 0;
        tAudTargets = round(double(msExptInfo(iexp).tAudTargets),2,'significant');
        audTrials = tAudTargets > 0;
        nTrials = length(tVisTargets);
        hit = msExptInfo(iexp).hit;
        miss = msExptInfo(iexp).miss;
        
        tInv = round(double(msExptInfo(iexp).tInvTargets),2,'significant');
        invTrials = tInv > 0;
        invHit = msExptInfo(iexp).invHit;
        invMiss = msExptInfo(iexp).invMiss;
        invTargets = unique(tInv);
        invTargets = invTargets(2:end);
        if strcmp(msExptInfo(iexp).invType,'vis')
            tInvVisTargets = tInv;
            tInvAudTargets = nan(1,nTrials);
            invHRexpt{iexp} = nan(1,length(invTargets));
            valHRexpt{iexp} = nan(1,length(invTargets));
            for i = 1:length(invTargets)
                invHRexpt{iexp}(i) = sum(invHit & tInv == invTargets(i)) ./ ...
                    sum((invHit | invMiss) & tInv == invTargets(i));                
                valHRexpt{iexp}(i) = sum(hit & tVisTargets == invTargets(i)) ./ ...
                    sum((hit | miss) & tVisTargets == invTargets(i));
            end
        elseif strcmp(msExptInfo(iexp).invType,'aud')
            tInvVisTargets = nan(1,nTrials);
            tInvAudTargets = tInv;
            invHRexpt{iexp} = nan(1,length(invTargets));
            valHRexpt{iexp} = nan(1,length(invTargets));
            for i = 1:length(invTargets)
                invHRexpt{iexp}(i) = sum(invHit & tInv == invTargets(i)) ./ ...
                    sum((invHit | invMiss) & tInv == invTargets(i));                
                valHRexpt{iexp}(i) = sum(hit & tAudTargets == invTargets(i)) ./ ...
                    sum((hit | miss) & tAudTargets == invTargets(i));
            end
        end

        
        trLength = msExptInfo(iexp).trLength;
        invTrLength = msExptInfo(iexp).invTrLength;
        
        fa = msExptInfo(iexp).fa;
        minCycMs = (msExptInfo(iexp).stimOnTime+msExptInfo(iexp).stimOffTime)...
            *minCyclesFA;
        trLengthMs = msExptInfo(iexp).trLengthMs;
        longEnoughTrials = trLengthMs > (minCycMs+100);
        nFACyc = trLength - minCyclesFA;
        FAR(visualTrials,iexp) = ...
            sum(fa(longEnoughTrials & visTrials & ~invTrials))/...
            sum(nFACyc(longEnoughTrials & visTrials & ~invTrials));
        FAR(auditoryTrials,iexp) = ...
            sum(fa(longEnoughTrials & audTrials & ~invTrials))/...
            sum(nFACyc(longEnoughTrials & audTrials & ~invTrials));
        

        if isfield(msCmlvData,'tVisTargets')
            msCmlvData.tVisTargets = cat(2,msCmlvData.tVisTargets,tVisTargets);
            msCmlvData.tAudTargets = cat(2,msCmlvData.tAudTargets,tAudTargets);
            msCmlvData.tInvVisTargets = cat(2,msCmlvData.tInvVisTargets,tInvVisTargets);
            msCmlvData.tInvAudTargets = cat(2,msCmlvData.tInvAudTargets,tInvAudTargets);
            msCmlvData.hit = cat(2,msCmlvData.hit,hit);
            msCmlvData.miss = cat(2,msCmlvData.miss,miss);
            msCmlvData.invHit = cat(2,msCmlvData.invHit,invHit);
            msCmlvData.invMiss = cat(2,msCmlvData.invMiss,invMiss);
            msCmlvData.trLength = cat(2,msCmlvData.trLength,trLength);
            msCmlvData.invTrLength = cat(2,msCmlvData.invTrLength,invTrLength);
        else
            msCmlvData.tVisTargets = tVisTargets;
            msCmlvData.tAudTargets = tAudTargets;
            msCmlvData.tInvVisTargets = tInvVisTargets;
            msCmlvData.tInvAudTargets = tInvAudTargets;
            msCmlvData.hit = hit;
            msCmlvData.miss = miss;
            msCmlvData.invHit = invHit;
            msCmlvData.invMiss = invMiss;
            msCmlvData.trLength = trLength;
            msCmlvData.invTrLength = invTrLength;
        end
    end

    %%
    HR_lim = [0 110];
    HR_label = [0:20:110];
    cueColors = {'k';'c'};
    cueNames = {'Valid';'Invalid'};
    fa_lim = [0 15];
    subHR_lim = [-110 110];
    subHR_label = -100:20:100;
    %%
    invType = {msExptInfo.invType};
    mFAR = mean(FAR,1)*100;
    
    repFAR = cell2mat(cellfun(@(x,y) repmat(x,1,length(y)),...
        num2cell(mFAR),valHRexpt,'unif',0));
    valHRAll = cell2mat(valHRexpt)*100;
    invHRAll = cell2mat(invHRexpt)*100;
    
    figure
    suptitle(mouseName)
    subplot 131
    scatter(repFAR,valHRAll)
    figXAxis([],'False Alarm Rate (%)',fa_lim)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Valid Trials')
    subplot 132
    scatter(repFAR,invHRAll)
    figXAxis([],'False Alarm Rate (%)',fa_lim)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Invalid Trials')
    subplot 133
    scatter(repFAR,valHRAll - invHRAll)
    figXAxis([],'False Alarm Rate (%)',fa_lim)
    figYAxis([],'Vaild - Invalid Hit Rate (%)',subHR_lim,subHR_label,subHR_label)
    figAxForm
    title('Invalid Trials')  
    
    print(fullfile(fnout,[fileName '_FAR']),'-dpdf','-fillpage')
    
    figure
    scatter(valHRAll,invHRAll,'ko')
    hold on
    plot(HR_label,HR_label,'k--')
    h = errorbar(mean(valHRAll),nanmean(invHRAll),ste(invHRAll,2),ste(invHRAll,2),...
        ste(valHRAll,2),ste(valHRAll,2),'ko');
    h.MarkerFaceColor = 'k'
    figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title(mouseName)
    
    print(fullfile(fnout,[fileName '_scatter_valInvHR']),'-dpdf','-fillpage')
    
    
    %%
    visTargets = unique(msCmlvData.tVisTargets);
    visTargets = visTargets(2:end);
    audTargets = unique(msCmlvData.tAudTargets);
    audTargets = audTargets(2:end);

    visLevels_lim = [min(visTargets)-1 110];
    visLevels_label = [10 100];
    audLevels_lim = [min(audTargets)-(0.5*min(audTargets)) 1.1];
    audLevels_label = [0 0.001 0.01 0.1 1];

    %%
    visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
    audBinEdges = exp(linspace(log(min(audTargets)-(0.5*min(audTargets))),...
        log(max(audTargets)),nBins+1));
    [~,~,visBinInd] = histcounts(msCmlvData.tVisTargets,visBinEdges);
    [~,~,audBinInd] = histcounts(msCmlvData.tAudTargets,audBinEdges);
    [~,~,invVisBinInd] = histcounts(msCmlvData.tInvVisTargets,visBinEdges);
    [~,~,invAudBinInd] = histcounts(msCmlvData.tInvAudTargets,audBinEdges);

    %%
    nVisHits = nan(1,nBins);
    nVisMisses = nan(1,nBins);
    nAudHits = nan(1,nBins);
    nAudMisses = nan(1,nBins);
    visTargetsBinned = nan(1,nBins);
    audTargetsBinned = nan(1,nBins);
    visTargetsSte = nan(1,nBins);
    audTargetsSte = nan(1,nBins);
    nInvVisHits = nan(1,nBins);
    nInvVisMisses = nan(1,nBins);
    nInvAudHits = nan(1,nBins);
    nInvAudMisses = nan(1,nBins);
    invVisTargetsBinned = nan(1,nBins);
    invAudTargetsBinned = nan(1,nBins);
    invVisTargetsSte = nan(1,nBins);
    invAudTargetsSte = nan(1,nBins);
    for ibin = 1:nBins
        ind = (msCmlvData.hit | msCmlvData.miss) & visBinInd == ibin;
        if sum(ind) > minTrN_ms
        nVisHits(ibin) = sum(msCmlvData.hit & visBinInd == ibin);
        nVisMisses(ibin) = sum(msCmlvData.miss & visBinInd == ibin);
        visTargetsBinned(ibin) = mean(msCmlvData.tVisTargets(ind));
        visTargetsSte(ibin) = ste(msCmlvData.tVisTargets(ind),2);
        end
        ind = (msCmlvData.hit | msCmlvData.miss) & audBinInd == ibin;
        if sum(ind) > minTrN_ms
        nAudHits(ibin) = sum(msCmlvData.hit & audBinInd == ibin);
        nAudMisses(ibin) = sum(msCmlvData.miss & audBinInd == ibin);
        audTargetsBinned(ibin) = mean(msCmlvData.tAudTargets(ind));
        audTargetsSte(ibin) = ste(msCmlvData.tAudTargets(ind),2);
        end
        ind = (msCmlvData.invHit | msCmlvData.invMiss) & invVisBinInd == ibin;
        if sum(ind) > minTrN_ms
            nInvVisHits(ibin) = sum(invVisBinInd == ibin & msCmlvData.invHit);
            nInvVisMisses(ibin) = sum(invVisBinInd == ibin & msCmlvData.invMiss);
            invVisTargetsBinned(ibin) = mean(msCmlvData.tInvVisTargets(ind));
            invVisTargetsSte(ibin) = ste(msCmlvData.tInvVisTargets(ind),2);
        end
        ind = (msCmlvData.invHit | msCmlvData.invMiss) & invAudBinInd == ibin;
        if sum(ind) > minTrN_ms  
            nInvAudHits(ibin) = sum(invAudBinInd == ibin & msCmlvData.invHit);
            nInvAudMisses(ibin) = sum(invAudBinInd == ibin & msCmlvData.invMiss);
            invAudTargetsBinned(ibin) = mean(msCmlvData.tInvAudTargets(ind));
            invAudTargetsSte(ibin) = ste(msCmlvData.tInvAudTargets(ind),2);
        end
    end
    visInd = ~isnan(nVisHits);
    invVisInd = ~isnan(nInvVisHits); 
    audInd = ~isnan(nAudHits);
    invAudInd = ~isnan(nInvAudHits); 

    [visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
    [audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
    [invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
        nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
    [invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
        nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));

    %% 
    nTrials = nVisHits(visInd)+nVisMisses(visInd);
    msVisFit = weibullFitLG(visTargetsBinned(visInd), visHR, 1,0, {'nTrials',nTrials});
    
    maxI = max(visTargetsBinned(visInd));
    minI = min(visTargetsBinned(visInd));
    msVisXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    
    nTrials = nAudHits(audInd)+nAudMisses(audInd);
    msAudFit = weibullFitLG(audTargetsBinned(audInd), audHR, 1,0, {'nTrials',nTrials});
    
    maxI = max(audTargetsBinned(audInd));
    minI = min(audTargetsBinned(audInd));
    msAudXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    
    %%
    bxSummaryTable = table(nExpAnalyzed{visualTrials},nansum(nVisHits+nVisMisses),...
        nansum(nInvVisHits+nInvVisMisses),nExpAnalyzed{auditoryTrials},...
        nansum(nAudHits+nAudMisses),nansum(nInvAudHits+nInvAudMisses));
    bxSummaryTable.Properties.VariableNames = {'nVisSessions','nVisValidIncluded',....
        'nVisInvalidIncl','nAudSessions','nAudValidInlcuded','nAudInvalidIncluded'};

    %%
    figure
    suptitle(mouseName)
    subplot 131
    x = visTargetsBinned(visInd);
    xErr = visTargetsSte(visInd);
    y = visHR.*100;
    yErrL = y - (visHR95ci(:,1)'.*100);
    yErrU = (visHR95ci(:,2)'.*100) - y;
    h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
    h.Color = cueColors{valid};
    h.MarkerFaceColor = [1 1 1];
    hold on
    x = invVisTargetsBinned(invVisInd);
    xErr = invVisTargetsSte(invVisInd);
    y = invVisHR.*100;
    yErrL = y - (invVisHR95ci(:,1)'.*100);
    yErrU = (invVisHR95ci(:,2)'.*100) - y;
    h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
    h.Color = cueColors{invalid};
    h.MarkerFaceColor = [1 1 1];
    hold on
    h=line(msVisXGrid, msVisFit.modelFun(msVisFit.coefEsts, msVisXGrid)*100);
    h.Color = 'k';
    vline(msVisFit.thresh,'k--',sprintf('threshold = %s',num2str(msVisFit.thresh)))
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Visual Trials')

    subplot 132
    x = audTargetsBinned(audInd);
    xErr = audTargetsSte(audInd);
    y = audHR.*100;
    yErrL = y - (audHR95ci(:,1)'.*100);
    yErrU = (audHR95ci(:,2)'.*100) - y;
    h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
    h.Color = cueColors{valid};
    h.MarkerFaceColor = [1 1 1];
    hold on
    x = invAudTargetsBinned(invAudInd);
    xErr = invAudTargetsSte(invAudInd);
    y = invAudHR.*100;
    yErrL = y - (invAudHR95ci(:,1)'.*100);
    yErrU = (invAudHR95ci(:,2)'.*100) - y;
    h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
    h.Color = cueColors{invalid};
    h.MarkerFaceColor = [1 1 1];
    hold on
    h=line(msAudXGrid, msAudFit.modelFun(msAudFit.coefEsts, msAudXGrid)*100);
    h.Color = 'k';
    vline(msAudFit.thresh,'k--',sprintf('threshold = %s',num2str(msAudFit.thresh)))
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Auditory Trials')

    subplot 133
    ax = gca;
    ax.Visible = 'off';
    h = text(0.6,0.5,bxSummaryTable.Properties.VariableNames',...
        'FontWeight', 'bold', 'FontSize', 10,'HorizontalAlignment','right');
    bxSummaryNs = cellfun(@num2str,...
        {nExpAnalyzed{visualTrials},nansum(nVisHits+nVisMisses),...
        nansum(nInvVisHits+nInvVisMisses),nExpAnalyzed{auditoryTrials},...
        nansum(nAudHits+nAudMisses),nansum(nInvAudHits+nInvAudMisses)},...
        'unif',0)';
    text(0.7, 0.5,bxSummaryNs,'FontWeight', 'bold', 'FontSize', 10);

    print(fullfile(fnout,fileName),'-dpdf','-fillpage')
    save(fullfile(fnout,[fileName '_data']),'msExptInfo')
end