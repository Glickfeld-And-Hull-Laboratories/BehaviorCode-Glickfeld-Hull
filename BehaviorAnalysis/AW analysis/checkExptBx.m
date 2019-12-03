clear all
close all

ds = 'FSV_V1';
eval(ds);
audControl = false;
exptInd = 38;
if ~isnan(exptInd)
    expt = expt(exptInd);
end

% nmice = length(ms2analyze);
bxParams_FSAV_attnV1ms

rc = behavConstsAV;
nexp = size(expt,2);
%%
bxStruct = createBxStruct(expt);

for iexp = 1:nexp
    smoothHRBin = 10;
    
    if audControl
        easyTrialInd = bxStruct(iexp).tVisTargets == max(bxStruct(iexp).tVisTargets);
    else
        easyTrialInd = bxStruct(iexp).tVisTargets == max(bxStruct(iexp).tVisTargets);
    end
    
    h = bxStruct(iexp).hit(easyTrialInd);
    m = bxStruct(iexp).miss(easyTrialInd);
    movHR = movsum(m,smoothHRBin)./(movsum(m,smoothHRBin)+movsum(h,smoothHRBin));
    
    e = bxStruct(iexp).fa;
    movEarly = movsum(e,smoothHRBin)./smoothHRBin;
    
    ntrials = length(e);
    
    figure
    subplot 121
    hold on
    plot(find(easyTrialInd),movHR,'m','LineWidth',2)
    plot(1:ntrials,movEarly,'c','LineWidth',2)
    figXAxis([],'Trial Number',[1 ntrials])
    figYAxis([],'Lapse (Pink) or Early (Cyan) Rate',[0 1])
    figAxForm
    if length(bxStruct(iexp).trialsPerRun) > 1
        runEndTrial = cumsum(bxStruct(iexp).trialsPerRun(1:end-1))+1;
        vline(runEndTrial,'k:')
    end
    title([bxStruct(iexp).sn '-' bxStruct(iexp).date])
    
    visTargets = unique(bxStruct(iexp).tVisTargets);
    if audControl
        audTargets = visTargets;
    else
        audTargets = unique(bxStruct(iexp).tAudTargets);
    end
    nvis = length(visTargets);
    naud = length(audTargets);
    
    nHits = nan(2,max([nvis,naud]));
    nMisses = nan(2,max([nvis,naud]));
    for i = 1:nvis
        nHits(visualTrials,i) = sum(bxStruct(iexp).tVisTargets == visTargets(i) & ...
            bxStruct(iexp).hit & bxStruct(iexp).isTrueVisTrial);
        nMisses(visualTrials,i) = sum(bxStruct(iexp).tVisTargets == visTargets(i) & ...
            bxStruct(iexp).miss & bxStruct(iexp).isTrueVisTrial);
    end
    for i = 1:naud
        if audControl
            nHits(auditoryTrials,i) = sum(bxStruct(iexp).tVisTargets == visTargets(i) & ...
                bxStruct(iexp).hit & ~bxStruct(iexp).isTrueVisTrial);
            nMisses(auditoryTrials,i) = sum(bxStruct(iexp).tVisTargets == visTargets(i) & ...
                bxStruct(iexp).miss & ~bxStruct(iexp).isTrueVisTrial);
        else
            nHits(auditoryTrials,i) = sum(bxStruct(iexp).tAudTargets == audTargets(i) & ...
                bxStruct(iexp).hit);
            nMisses(auditoryTrials,i) = sum(bxStruct(iexp).tAudTargets == audTargets(i) & ...
                bxStruct(iexp).miss);
        end
    end
    ind = ~isnan(nHits(visualTrials,:));
    [visHR,visCI] = binofit(nHits(visualTrials,ind),...
        nHits(visualTrials,ind)+nMisses(visualTrials,ind));    
    visInd = ~isnan(visHR);
    ind = ~isnan(nHits(auditoryTrials,:));
    [audHR,audCI] = binofit(nHits(auditoryTrials,ind),...
        nHits(auditoryTrials,ind)+nMisses(auditoryTrials,ind));
    audInd = ~isnan(audHR);
    
    subplot 122
    hold on
%     h = errorbar(visTargets(visInd),visHR(visInd),...
%         visHR(visInd)-visCI(visInd,1)',visCI(visInd,2)'-visHR(visInd),...
%         '.-','MarkerSize',10,'LineWidth',2);
    h = plot(visTargets(visInd),visHR(visInd),...
        '.-','MarkerSize',10,'LineWidth',2);
    h.Color = cueColor{visualTrials};
%     h = errorbar(audTargets(audInd),audHR(audInd),...
%         audHR(audInd)-audCI(audInd,1)',audCI(audInd,2)'-audHR(audInd),...
%         '.-','MarkerSize',10,'LineWidth',2);
     h = plot(audTargets(audInd).*100,audHR(audInd),...
        '.-','MarkerSize',10,'LineWidth',2);
    h.Color = cueColor{auditoryTrials};
    figXAxis([],'Orientation',[0 100])
    figYAxis([],'Hit Rate',[0 1])
    figAxForm
    hline(0.9,'k:')
    title(sprintf('Early Rate: %s',num2str(round(bxStruct(iexp).pctEarly,2,'significant'))))
end