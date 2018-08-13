% run eaMsBxSummary and check bxParams_FSAV before starting

%% load data
clear all
close all
rc = behavConstsAV;
exptSummaryDir = fullfile(rc.ashley,...
    'Manuscripts','Attention V1','Mouse Info.xlsx');
exptSummaryInfo = readtable(exptSummaryDir);
fnout = fullfile(rc.ashley,'Manuscripts','Attention V1','Matlab Figs');

ms2analyze = cellfun(@num2str,num2cell(exptSummaryInfo.SubjectNumber),...
    'unif',0);
% ms2analyze = exptSummaryInfo.SubjectNumber';
nMice = length(ms2analyze);

exampleMouse = '668';
exMsInd = find(strcmp(ms2analyze,exampleMouse));

bxParams_FSAV
%% compile data
msSumStruct = [];
msExptStruct = cell(1,nMice);
for im = 1:nMice
    mouseName = ms2analyze{im};
    fn = fullfile(rc.ashleyAnalysis,mouseName,'behavior');
    load(fullfile(fn,[mouseName,'bxSummary_dataAnalyzed']))
    msSumStruct = cat(1,msSumStruct,msCmlvData);
    if im == exMsInd
        exMsExptInfo = msExptAnalyzed;
    end
end

msHR = struct;
HR = cell(2,2);
HR(:) = {nan(nMice,nBins)};
targets = cell(2,2);
targets(:) = {nan(nMice,nBins)};
for im = 1:nMice
    
    msCmlvData = msSumStruct(im);
    
    visTargets = unique(msCmlvData.tVisTargets);
    visTargets = visTargets(2:end);
    audTargets = unique(msCmlvData.tAudTargets);
    audTargets = audTargets(2:end);

    visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
    audBinEdges = exp(linspace(...
        log(min(audTargets(audTargets > 0.00001))-...
        (0.5*min(audTargets(audTargets > 0.00001)))),...
        log(max(audTargets)),nBins+1));
    [~,~,visBinInd] = histcounts(msCmlvData.tVisTargets,visBinEdges);
    [~,~,audBinInd] = histcounts(msCmlvData.tAudTargets,audBinEdges);
    [~,~,invVisBinInd] = histcounts(msCmlvData.tInvVisTargets,visBinEdges);
    [~,~,invAudBinInd] = histcounts(msCmlvData.tInvAudTargets,audBinEdges);
    
    nVisInv(im) = sum(msCmlvData.tInvVisTargets > 0);
    nAudInv(im) = sum(msCmlvData.tInvAudTargets > 0);
        
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
    
    [visAttnTest,visHRall,visInvHRall] = allTrialsAttnTest(msCmlvData.tVisTargets,...
        msCmlvData.tInvVisTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);
    [audAttnTest,audHRall,audInvHRall] = allTrialsAttnTest(msCmlvData.tAudTargets,...
        msCmlvData.tInvAudTargets,msCmlvData.hit,msCmlvData.miss,...
        msCmlvData.invHit,msCmlvData.invMiss);

    [visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
    [audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
    if sum(nInvVisHits(invVisInd)+nInvVisMisses(invVisInd)) > 0
        [invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
            nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
    else
        invVisHR = nan;
        invVisHR95ci = nan(1,2);
    end
    if sum(nInvAudHits(invAudInd)+nInvAudMisses(invAudInd)) > 0
        [invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
            nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));
    else
        invAudHR = nan;
        invAudHR95ci = nan(1,2);
    end
    
    nTrials = nVisHits(visInd)+nVisMisses(visInd);
    msVisFit = weibullFitLG(visTargetsBinned(visInd), visHR, 0,0, {'nTrials',nTrials});

    maxI = max(visTargetsBinned(visInd));
    minI = min(visTargetsBinned(visInd));
    msVisXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);

    nTrials = nAudHits(audInd)+nAudMisses(audInd);
    msAudFit = weibullFitLG(audTargetsBinned(audInd), audHR, 0,0, {'nTrials',nTrials});

    maxI = max(audTargetsBinned(audInd));
    minI = min(audTargetsBinned(audInd));
    msAudXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
    
    
%     f = msVisFit.modelFun(msVisFit.coefEsts, msVisXGrid);
%     oriAtHighThresh = msVisXGrid(find(f > highThreshold,1));
%     f = msAudFit.modelFun(msAudFit.coefEsts, msAudXGrid);
%     ampAtHighThresh = msAudXGrid(find(f > highThreshold,1));
%     lowVisAttnTest = belowThreshAttnTest(oriAtHighThresh,msCmlvData.tVisTargets,...
%         msCmlvData.tInvVisTargets,msCmlvData.hit,msCmlvData.miss,...
%         msCmlvData.invHit,msCmlvData.invMiss);
%     lowAudAttnTest = belowThreshAttnTest(ampAtHighThresh,msCmlvData.tAudTargets,...
%         msCmlvData.tInvAudTargets,msCmlvData.hit,msCmlvData.miss,...
%         msCmlvData.invHit,msCmlvData.invMiss);
        
    [valHR_highThreshold_vis, invHR_highThreshold_vis] = ...
        getMatchedHighThresholdHR(visTargetsBinned,msVisFit,highThreshold,...
        msCmlvData.tVisTargets,msCmlvData.tInvVisTargets,...
        msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
    
    [valHR_highThreshold_aud, invHR_highThreshold_aud] = ...
        getMatchedHighThresholdHR(audTargetsBinned,msAudFit,highThreshold,...
        msCmlvData.tAudTargets,msCmlvData.tInvAudTargets,...
        msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
    
    msHR(im).av(visualTrials).cue(valid).HR = visHR.*100;
    msHR(im).av(visualTrials).cue(valid).HR95CI = visHR95ci.*100;
    msHR(im).av(visualTrials).cue(valid).targets = visTargetsBinned(visInd);
    msHR(im).av(visualTrials).cue(valid).targetsErr = visTargetsSte(visInd);
    msHR(im).av(visualTrials).cue(valid).fit = msVisFit;
    msHR(im).av(visualTrials).cue(valid).fitGrid = msVisXGrid;
    msHR(im).av(visualTrials).cue(valid).hiLoHR = valHR_highThreshold_vis;
    msHR(im).av(visualTrials).attnTest = visAttnTest;   
    msHR(im).av(visualTrials).matchedHRall = [visHRall,visInvHRall];
%     msHR(im).av(visualTrials).threshAttnTest = lowVisAttnTest;    
    
    HR{visualTrials,valid}(im,visInd) = ...
        msHR(im).av(visualTrials).cue(valid).HR;
    targets{visualTrials,valid}(im,visInd) = ...
        msHR(im).av(visualTrials).cue(valid).targets;
    
    msHR(im).av(auditoryTrials).cue(valid).HR = audHR.*100;
    msHR(im).av(auditoryTrials).cue(valid).HR95CI = audHR95ci.*100;
    msHR(im).av(auditoryTrials).cue(valid).targets = audTargetsBinned(audInd);
    msHR(im).av(auditoryTrials).cue(valid).targetsErr = audTargetsSte(audInd);
    msHR(im).av(auditoryTrials).cue(valid).fit = msAudFit;
    msHR(im).av(auditoryTrials).cue(valid).fitGrid = msAudXGrid;
    msHR(im).av(auditoryTrials).cue(valid).hiLoHR = valHR_highThreshold_aud;
    msHR(im).av(auditoryTrials).attnTest = audAttnTest;
    msHR(im).av(auditoryTrials).matchedHRall = [audHRall,audInvHRall];
%     msHR(im).av(auditoryTrials).threshAttnTest = lowAudAttnTest;
    
    
    HR{auditoryTrials,valid}(im,audInd) = ...
        msHR(im).av(auditoryTrials).cue(valid).HR;
    targets{auditoryTrials,valid}(im,audInd) = ...
        msHR(im).av(auditoryTrials).cue(valid).targets;
    
    msHR(im).av(visualTrials).cue(invalid).HR = invVisHR.*100;
    msHR(im).av(visualTrials).cue(invalid).HR95CI = invVisHR95ci.*100;
    msHR(im).av(visualTrials).cue(invalid).targets = invVisTargetsBinned(invVisInd);
    msHR(im).av(visualTrials).cue(invalid).targetsErr = invVisTargetsSte(invVisInd);
    msHR(im).av(visualTrials).cue(invalid).hiLoHR = invHR_highThreshold_vis;
    
    HR{visualTrials,invalid}(im,invVisInd) = ...
        msHR(im).av(visualTrials).cue(invalid).HR;
    targets{visualTrials,invalid}(im,invVisInd) = ...
        msHR(im).av(visualTrials).cue(invalid).targets;
    
    msHR(im).av(auditoryTrials).cue(invalid).HR = invAudHR.*100;
    msHR(im).av(auditoryTrials).cue(invalid).HR95CI = invAudHR95ci.*100;
    msHR(im).av(auditoryTrials).cue(invalid).targets = invAudTargetsBinned(invAudInd);
    msHR(im).av(auditoryTrials).cue(invalid).targetsErr = invAudTargetsSte(invAudInd);
    msHR(im).av(auditoryTrials).cue(invalid).hiLoHR = invHR_highThreshold_aud;
    
    HR{auditoryTrials,invalid}(im,invAudInd) = ...
        msHR(im).av(auditoryTrials).cue(invalid).HR;
    targets{auditoryTrials,invalid}(im,invAudInd) = ...
        msHR(im).av(auditoryTrials).cue(invalid).targets;    
end

%%
visAttnTest = nan(1,nMice);
audAttnTest = nan(1,nMice);
for i = 1:nMice
    visAttnTest(i) = msHR(i).av(visualTrials).attnTest;
    audAttnTest(i) = msHR(i).av(auditoryTrials).attnTest;
end

%%

HR_lim = [0 110];
HR_label = [0:20:110];

visLevels_lim = [min(targets{visualTrials,valid}(:))-1 110];
visLevels_label = [11.25 22.5 45 90];
audLevels_lim = [min(targets{auditoryTrials,valid}(:))-...
    (0.5*min(audTargets(audTargets > 0.00001))) 1.1];
audLevels_label = [0 0.001 0.01 0.1 1];

%% plot HR summary
setFigParams4Print('portrait')
set(0,'defaultAxesFontSize',10)
fitsOris = cell(1,2);
fitsOris(:) = {nan(100,nMice)};
fitsHR = cell(1,2);
fitsHR(:) = {nan(100,nMice)};
HRfig = figure;
for im = 1:nMice
    for iav = 1:2
        for icue = 1:2
            x = msHR(im).av(iav).cue(icue).targets;
            xerr = msHR(im).av(iav).cue(icue).targetsErr;
            y = msHR(im).av(iav).cue(icue).HR;
            ylerr = y - (msHR(im).av(iav).cue(icue).HR95CI(:,1)');
            yuerr = (msHR(im).av(iav).cue(icue).HR95CI(:,2)') - y;
            if icue == valid
                f = msHR(im).av(iav).cue(icue).fit;
                fitX = msHR(im).av(iav).cue(icue).fitGrid;
                fitY = f.modelFun(f.coefEsts, fitX).*100;
                fitsOri{iav}(:,im) = fitX;
                fitsHR{iav}(:,im) = fitY;
               subplot(5,2,iav+2)
                hold on
                h = plot(fitX,fitY,'-');
                h.Color = cueColor{icue};
            elseif ~all(isnan(y))
               subplot(5,2,iav+2)
                hold on
                h = plot(x,y,'.-');
                h.Color = cueColor{icue};
            end 
            if ~all(isnan(y))
               subplot(5,2,iav)
                hold on
                h = plot(x,y,'.-');
                h.Color = cueColor{icue};
            end
        end
    end
end 

figure(HRfig);
for iav = 1:2
    for icue = 1:2
        x = nanmean(targets{iav,icue},1);
        xerr = ste(targets{iav,icue},1);
        y = nanmean(HR{iav,icue},1);
        yerr = ste(HR{iav,icue},1);
        ind = ~isnan(y);
        
       subplot(5,2,iav)
        hold on
        h = errorbar(x(ind),y(ind),yerr(ind),yerr(ind),xerr(ind),xerr(ind),'o-');
        h.Color = cueColor{icue};
        h.LineWidth = 1;
        h.MarkerFaceColor = [1 1 1];
        
        if icue == valid
            fitX = mean(fitsOri{iav},2);
            fitY = mean(fitsHR{iav},2);
            fitYerr = ste(fitsHR{iav},2);
           subplot(5,2,iav+2)
            hold on
            h = shadedErrorBar(fitX,fitY,fitYerr,'-');
            h.mainLine.Color = cueColor{icue};
            h.mainLine.LineWidth = 2;
            h.patch.FaceColor = cueColor{icue} +((1 - cueColor{icue}).*0.5);
            h.edge(1).Color = h.patch.FaceColor;
        else
           subplot(5,2,iav+2)
            hold on
%             h = errorbar(x(ind),y(ind),yerr(ind),yerr(ind),xerr(ind),xerr(ind),'o-');
            h = shadedErrorBar(x(ind),y(ind),yerr(ind));
            h.mainLine.Color = cueColor{icue};
            h.mainLine.LineWidth = 2;
            h.patch.FaceColor = cueColor{icue} +((1 - cueColor{icue}).*0.5);
            h.edge(1).Color = h.patch.FaceColor;
        end
    end
end

hiLoHR = cell(2,2);
hiLoHR(:) = {nan(nMice,2)};
allHR = cell(1,2);
allHR(:) = {nan(nMice,2)};
for im = 1:nMice
    for iav = 1:2
        for ithresh = 1:2
            x = 1:2;
            y = [msHR(im).av(iav).cue(1).hiLoHR(ithresh),...
                msHR(im).av(iav).cue(2).hiLoHR(ithresh)].*100;
            hiLoHR{iav,ithresh}(im,:) = y;            
            if ithresh == 1
                subplot(5,2,iav+4)
                hold on
                title(sprintf('Targets < %s Threshold',num2str(highThreshold.*100)))
            elseif ithresh == 2
                subplot(5,2,iav+6)
                hold on
                title(sprintf('Targets > %s Threshold',num2str(highThreshold.*100)))
            end
            h = plot(x,y,'-');
            h.Color = hiLoColor{ithresh};
            if im == nMice
                y = mean(hiLoHR{iav,ithresh},1);
                yerr = ste(hiLoHR{iav,ithresh},1);
                h = errorbar(x,y,yerr,'o-');
                h.LineWidth = 1;
                h.Color = hiLoColor{ithresh};
                h.MarkerFaceColor = [1 1 1];
                figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
                figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
                figAxForm
            end
            
        end
        
        x = 1:2;
        y = msHR(im).av(iav).matchedHRall.*100;
        p = msHR(im).av(iav).attnTest;
        allHR{iav}(im,:) = y;
        subplot(5,2,iav+8)
        hold on
        h = plot(x,y,'-');
        if p < attnTestAlpha
            h.Color = 'k';
        else
            h.Color = [0.5 0.5 0.5];
        end
        if im == nMice
            y = mean(allHR{iav},1);
            yerr = ste(allHR{iav},1);
            h = errorbar(x,y,yerr,'o-');
            h.LineWidth = 1;
            h.Color = 'k';
            h.MarkerFaceColor = [1 1 1];
            figXAxis([],'',[0 3],x,{'Vaild';'Invalid'})
            figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
            figAxForm
            title(sprintf('All Trials, alpha = %s',num2str(round(attnTestAlpha,2,'significant'))))
        end
    end
end


for iplot = 1:2
    plotOffset = (iplot-1)*2;
   subplot(5,2,visualTrials+plotOffset)
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Visual Trials')
   hold on
   hline(highThreshold*100,'k:')

   subplot(5,2,auditoryTrials+plotOffset)
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Auditory Trials')
   hold on
   hline(highThreshold*100,'k:')
end
print(fullfile(fnout,'bxSummary_allMice'),'-dpdf','-fillpage')
%% example mouse

nexp = size(exMsExptInfo,2);
exExptHR = cell(2,2);
exExptHR(:) = {nan(2,nexp)};
exExptHighThresh = cell(2,2);
exExptHighThresh(:) = {nan(2,nexp)};
for i = 1:nexp 
    for iav = 1:2
        for icue = 1:2
            exExptHR{iav,icue}(:,i) = exMsExptInfo(i).av(iav).cue(icue).highThreshHR;
        end
    end
end

%% example mouse plot
exMsFig = figure;
suptitle(ms2analyze(exMsInd))
for iav = 1:2
    for icue = 1:2
        x = msHR(exMsInd).av(iav).cue(icue).targets;
        xerr = msHR(exMsInd).av(iav).cue(icue).targetsErr;
        y = msHR(exMsInd).av(iav).cue(icue).HR;
        ylerr = y -  msHR(exMsInd).av(iav).cue(icue).HR95CI(:,1)';
        yuerr = msHR(exMsInd).av(iav).cue(icue).HR95CI(:,2)' - y;
        
        subplot(2,2,iav)
        hold on
        h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'o-');
        h.Color = cueColor{icue};
        h.LineWidth = 1;
        h.MarkerFaceColor = [1 1 1];
        
        if icue == valid
            f = msHR(exMsInd).av(iav).cue(icue).fit;
            x = msHR(exMsInd).av(iav).cue(icue).fitGrid;
            y = f.modelFun(f.coefEsts, x).*100;
            subplot(2,2,iav)
            hold on
            h = plot(x,y,'-');
            h.LineWidth = 1;
            h.Color = cueColor{icue};
        end
    end
end
subplot(2,2,visualTrials)
ax = gca;
ax.XScale = 'log';
figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Visual Trials - Combined Expt')
   hold on
   hline(highThreshold*100,'k:')

subplot(2,2,auditoryTrials)
ax = gca;
ax.XScale = 'log';
figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Auditory Trials - Combined Expt')
   hold on
   hline(highThreshold*100,'k:')

for iav = 1:2
    for ithresh = 1:2
        subplot(2,2,iav+2)
        hold on
        x = exExptHR{iav,1}(ithresh,:).*100;
        y = exExptHR{iav,2}(ithresh,:).*100;
        h = plot(x,y,'o');
        h.Color = hiLoColor{ithresh};
        h.MarkerFaceColor = hiLoColor{ithresh};
        h = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'o');
        h.Color = hiLoColor{ithresh};
        h.MarkerFaceColor = [1 1 1];
        h.LineWidth = 1;
    end
    plot(HR_lim,HR_lim,'k--')
    figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    if iav == 1
        title('Visual Trials - Each Expt')
    elseif iav == 2
        title('Auditory Trials - Each Expt')
    end
end
print(fullfile(fnout,'bxSummary_exampleMouse'),'-dpdf','-fillpage')
%% 