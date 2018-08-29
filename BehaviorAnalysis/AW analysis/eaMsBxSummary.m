clear all
close all

ms2analyze = {'613';'614';'625';'668';'750'};
exampleDay = 10;
doExampleDay = false;
doLoadPreviousDataset = false;
doCheck4NewDates = true;
doRewarded = true;
doPlot = true;
%%
nmice = length(ms2analyze);
bxParams_FSAV

rc = behavConstsAV;

%%

HR_lim = [0 110];
HR_label = [0:20:110];
cueColors = {'k';'c'};
cueNames = {'Valid';'Invalid'};
threshColors = {'b';'k'};
threshNames = {'< Thresh';'> Thresh'};
fa_lim = [0 15];
subHR_lim = [-110 110];
subHR_label = -100:20:100;
react_lim = [0 600];
visThresh_lim = [10 60];
audThresh_lim = [0 0.2];
rt_lim = [100 500];
timeLim = [timeBins(1) - 100 timeBins(end)];
set(0,'defaultAxesFontSize',16)
%%
for im = 1:nmice
    mouseName = ms2analyze{im};
    fnout = fullfile(rc.ashleyAnalysis,mouseName,'behavior');
    if doRewarded
        fileName = sprintf('%s_bxSummary_%s',mouseName);
        savedDataName = [mouseName 'bxSummary_data'];
        xlsFileDir = rc.eaMouseIndexFilename;
        disp(savedDataName)
    else
        fileName = sprintf('%s_bxSummary_noRew_',mouseName);
        savedDataName = [mouseName 'bxSummary_noRew_data'];
        xlsFileDir = fullfile(rc.ashleyAnalysis,'Behavior',...
            'experimentIndexes','FSAV_noCatchRewEaMouseData.xls');
        disp(savedDataName)
    end
    if doLoadPreviousDataset
        load(fullfile(fnout,savedDataName))
        if doCheck4NewDates
            bxDataInfo = table2struct(readtable(xlsFileDir,'Sheet',mouseName));
            newDates = ~ismember({bxDataInfo.DateStr},{msExptInfo.date});
            if any(newDates)
                msExptInfo_newDates = bxFrmXLS(mouseName,bxDataInfo(newDates));
                msExptInfo = cat(2,msExptInfo,msExptInfo_newDates);
                save(fullfile(fnout,[mouseName 'bxSummary_data']),'msExptInfo')
            end
            if any(~ismember({msExptInfo.date},{bxDataInfo.DateStr}))
                dateInd = ismember({msExptInfo.date},{bxDataInfo.DateStr});
                msExptInfo = msExptInfo(dateInd);
                save(fullfile(fnout,savedDataName),'msExptInfo')
            end
        end
    else
        bxDataInfo = table2struct(readtable(xlsFileDir,'Sheet',mouseName));
        % fnin = rc.behavData;
        if ~exist(fnout,'dir')
            mkdir(fnout)
        end
        msExptInfo = bxFrmXLS(mouseName,bxDataInfo);
        save(fullfile(fnout,savedDataName),'msExptInfo')
    end
    %%
    disp(mouseName)
    nexp = size(msExptInfo,2);
    msCmlvData = struct;
    msExptAnalyzed = struct;    
    exptInd = false(1,nexp);
    for iexp = 1:nexp
        if iexp == 1
            exptN = 0;
        end
        if length(msExptInfo(iexp).visHR) > 1
            visHRCutoffPass = msExptInfo(iexp).visHR(end-1:end) < lapse_cutoff;
        else
            visHRCutoffPass = msExptInfo(iexp).visHR < lapse_cutoff;
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
        exptN = exptN+1;


        msExptAnalyzed(exptN).date = str2double(msExptInfo(iexp).date);
        tVisTargets = round(double(msExptInfo(iexp).tVisTargets),2,'significant');
        visTrials = tVisTargets > 0;
        tAudTargets = round(double(msExptInfo(iexp).tAudTargets),2,'significant');
        audTrials = tAudTargets > 0;
        hit = msExptInfo(iexp).hit;
        miss = msExptInfo(iexp).miss;
        valRT = msExptInfo(iexp).valReact;
        
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).nTrials = sum(tVisTargets > 0);
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).nTrials = sum(tAudTargets > 0);
        
        tInvVisTargets = round(double(msExptInfo(iexp).tInvVisTargets),2,'significant');
        tInvAudTargets = round(double(msExptInfo(iexp).tInvAudTargets),2,'significant');
        invTrials = tInvVisTargets > 0 | tInvAudTargets > 0;
        invHit = msExptInfo(iexp).invHit;
        invMiss = msExptInfo(iexp).invMiss;
        invRT = msExptInfo(iexp).invReact;
        
        cycLengthMs = double(msExptInfo(iexp).stimOnTime+msExptInfo(iexp).stimOffTime);
        tCyc = floor(msExptInfo(iexp).trLengthMs./cycLengthMs);
        lastStimRT = msExptInfo(iexp).trLengthMs - (tCyc.*cycLengthMs);
        tCyc = tCyc - 1;
        tCyc(tCyc < 0) = nan;
        oneStimBackRT = msExptInfo(iexp).trLengthMs - (tCyc.*cycLengthMs);
        tCyc = tCyc - 1;
        tCyc(tCyc < 0) = nan;
        twoStimBackRT = msExptInfo(iexp).trLengthMs - (tCyc.*cycLengthMs);
        
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).nTrials = sum(tInvVisTargets > 0);
        msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).nTrials = sum(tInvAudTargets > 0);
        
        invTargets = unique(tInvVisTargets(tInvVisTargets > 0));
        invVisHRexpt = nan(1,length(invTargets));
        valVisHRexpt = nan(1,length(invTargets)); 
        nInvExpt = nan(1,length(invTargets));
        nTrialsExpt = nan(1,length(invTargets));
        valVisRT = nan(1,length(invTargets));
        invVisRT = nan(1,length(invTargets));
        for i = 1:length(invTargets)
            invVisHRexpt(i) = sum(invHit & tInvVisTargets == invTargets(i)) ./ ...
                sum((invHit | invMiss) & tInvVisTargets == invTargets(i));                
            valVisHRexpt(i) = sum(hit & tVisTargets == invTargets(i)) ./ ...
                sum((hit | miss) & tVisTargets == invTargets(i));
            nInvExpt(i) = sum((invHit | invMiss) & tInvVisTargets == invTargets(i));
            nTrialsExpt(i) = ...
                sum((invHit | invMiss) & tInvVisTargets == invTargets(i)) + ...
                sum((hit | miss) & tVisTargets == invTargets(i));
            valVisRT(i) = mean(valRT(hit & tVisTargets == invTargets(i)));
            invVisRT(i) = mean(invRT(invHit & tInvVisTargets == invTargets(i)));
        end
        
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).targets = invTargets;
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).matchedHR = invVisHRexpt;
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).nMatchedTrials = nInvExpt;
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).matchedRT = invVisRT;
        msExptAnalyzed(exptN).av(visualTrials).nMatchedTrials = nTrialsExpt;
        
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).matchedHR = valVisHRexpt;
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).nMatchedTrials = nTrialsExpt - nInvExpt;
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).matchedRT = valVisRT;
        
        invTargets = unique(tInvAudTargets(tInvAudTargets > 0));
        invAudHRexpt = nan(1,length(invTargets));
        valAudHRexpt = nan(1,length(invTargets));
        nInvExpt = nan(1,length(invTargets));
        nTrialsExpt = nan(1,length(invTargets));
        valAudRT = nan(1,length(invTargets));
        invAudRT = nan(1,length(invTargets));
        for i = 1:length(invTargets)
            invAudHRexpt(i) = sum(invHit & tInvAudTargets == invTargets(i)) ./ ...
                sum((invHit | invMiss) & tInvAudTargets == invTargets(i));                
            valAudHRexpt(i) = sum(hit & tAudTargets == invTargets(i)) ./ ...
                sum((hit | miss) & tAudTargets == invTargets(i));
            nInvExpt(i) = ...
                sum((invHit | invMiss) & tInvAudTargets == invTargets(i));
            nTrialsExpt(i) = ...
                sum((invHit | invMiss) & tInvAudTargets == invTargets(i)) + ...
                sum((hit | miss) & tAudTargets == invTargets(i));
            valAudRT(i) = mean(valRT(hit & tAudTargets == invTargets(i)));
            invAudRT(i) = mean(invRT(invHit & tInvAudTargets == invTargets(i)));
        end

        msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).targets = invTargets;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).matchedHR = invAudHRexpt;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).nMatchedTrials = nInvExpt;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).matchedRT = invAudRT;
        msExptAnalyzed(exptN).av(auditoryTrials).nMatchedTrials = nTrialsExpt;
        
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).matchedHR = valAudHRexpt;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).nMatchedTrials = nTrialsExpt - nInvExpt;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).matchedRT = invAudRT;

        visStim = unique(tVisTargets(tVisTargets > 0));
        audStim = unique(tAudTargets(tAudTargets > 0));
        
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).targets = visStim;
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).targets = audStim;
        
        
        nVis = zeros(1,length(visStim));
        nAud = zeros(1,length(audStim));
        for i = 1:length(visStim)
            nVis(i) = sum((hit | miss) & tVisTargets == visStim(i));
        end
        for i = 1:length(audStim)
            nAud(i) = sum((hit | miss) & tAudTargets == audStim(i));
        end
        
        fa = msExptInfo(iexp).fa;
        stimOn = msExptInfo(iexp).stimOnTime;
        stimOff = msExptInfo(iexp).stimOffTime;
        minCycMs = (stimOn+stimOff).*minCyclesFA;
        trLengthMs = msExptInfo(iexp).trLengthMs;
        longEnoughTrials = trLengthMs > (minCycMs+100);
        tCycles = floor(trLengthMs./double(stimOn+stimOff));
        lastStimReactTime = trLengthMs - (tCycles.*double(stimOn+stimOff));
        trueFAInd = lastStimReactTime > 100 &  fa & longEnoughTrials;
        nFACyc = tCycles - minCyclesFA;
        
        nCycles = msExptInfo(iexp).nCycles;
        invCycles = msExptInfo(iexp).catchCyc;
        
        invTargetTimeMs = invCycles .* double(stimOn+stimOff);
        invTargetTimeMs(fa & ~invHit & ~invMiss) = nan;
        
        targetTimeMs = nCycles .* double(stimOn+stimOff);
        targetTimeMs(fa) = nan;
       
        msExptAnalyzed(exptN).av(visualTrials).falseAlarmRate = ...
            sum(fa(trueFAInd & visTrials & ~invTrials))/...
            sum(nFACyc(longEnoughTrials & visTrials & ~invTrials));
        nVisFATrials = sum(nFACyc(longEnoughTrials & visTrials & ~invTrials));
        msExptAnalyzed(exptN).av(auditoryTrials).falseAlarmRate = ...
            sum(fa(trueFAInd & audTrials & ~invTrials))/...
            sum(nFACyc(longEnoughTrials & audTrials & ~invTrials));
        nAudFATrials = sum(nFACyc(longEnoughTrials & audTrials & ~invTrials));
        
        HR = nan(1,length(visStim));
        for i = 1:length(visStim)
            HR(i) = sum(tVisTargets == visStim(i) & hit)./...
                sum((hit | miss) & tVisTargets == visStim(i));
        end
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).HR = HR;
        
        HR = nan(1,length(audStim));
        for i = 1:length(audStim)
            HR(i) = sum(tAudTargets == audStim(i) & hit)./...
                sum((hit | miss) & tAudTargets == audStim(i));
        end
        msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).HR = HR;
                
        visHR = msExptInfo(iexp).visHR;        
        ind = nVis > minTrN_expt;
%         msVisFit = weibullFitLG(visStim(ind),visHR(ind),1,0,{'nTrials',nVis(ind)});
        msVisFit = weibullFitLG([0 visStim(ind)],...
            [msExptAnalyzed(exptN).av(visualTrials).falseAlarmRate visHR(ind)]...
            ,1,0,{'nTrials',[nVisFATrials nVis(ind)]});

        msExptAnalyzed(exptN).av(visualTrials).threshold = msVisFit.thresh;
        
        [visHR_highThreshold, invVisHR_highThreshold, highThreshExpt,visMatchHighThreshTrialInd] = ...
            getMatchedHighThresholdHR(...
            visStim,msVisFit,highThreshold,tVisTargets,tInvVisTargets,...
            hit,miss,invHit,invMiss);
        
        msExptAnalyzed(exptN).av(visualTrials).highThreshold = highThreshExpt;
        msExptAnalyzed(exptN).av(visualTrials).cue(valid).highThreshHR = visHR_highThreshold;
        msExptAnalyzed(exptN).av(visualTrials).cue(invalid).highThreshHR = invVisHR_highThreshold;
        invVisHighThreshID = tInvVisTargets;
        invVisHighThreshID(tInvVisTargets > highThreshExpt) = 2;
        invVisHighThreshID(tInvVisTargets < highThreshExpt & tInvVisTargets > 0) = 1;
        visMatchedHighThreshID = visMatchHighThreshTrialInd(visMatchHighThreshTrialInd > 0 & (hit | miss));
        visMatchedHighThreshOutcome = hit(visMatchHighThreshTrialInd > 0 & (hit | miss));
        
        
        audHR = msExptInfo(iexp).audHR;
        ind = nAud > minTrN_expt;
        ind = ~isnan(audHR) & ind;
        if sum(ind) >= 2
            
%             msAudFit = weibullFitLG(audStim(ind),audHR(ind),1,0,{'nTrials',nAud(ind)});

            msAudFit = weibullFitLG([0 audStim(ind)],...
                [msExptAnalyzed(exptN).av(auditoryTrials).falseAlarmRate audHR(ind)],...
                1,0,{'nTrials',[nAudFATrials nAud(ind)]});        
            [audHR_highThreshold, invAudHR_highThreshold, highThreshExpt ,audMatchHighThreshTrialInd] = ...
                getMatchedHighThresholdHR(...
                audStim,msAudFit,highThreshold,tAudTargets,tInvAudTargets,...
                hit,miss,invHit,invMiss);
            audMatchedHighThreshID = audMatchHighThreshTrialInd(audMatchHighThreshTrialInd > 0 & (hit | miss));
            audMatchedHighThreshOutcome = hit(audMatchHighThreshTrialInd > 0 & (hit | miss)); 
            invAudHighThreshID = tInvAudTargets;
            invAudHighThreshID(tInvAudTargets > highThreshExpt) = 2;
            invAudHighThreshID(tInvAudTargets < highThreshExpt & tInvAudTargets > 0) = 1;  
            
            msExptAnalyzed(exptN).av(auditoryTrials).threshold = msAudFit.thresh;
            msExptAnalyzed(exptN).av(auditoryTrials).highThreshold = highThreshExpt;
            msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).highThreshHR = audHR_highThreshold;
            msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).highThreshHR = invAudHR_highThreshold;
        else
            audMatchedHighThreshID = [];
            audMatchedHighThreshOutcome = [];     
            invAudHighThreshID = [];
                        
            msExptAnalyzed(exptN).av(auditoryTrials).threshold = nan;
            msExptAnalyzed(exptN).av(auditoryTrials).highThreshold = nan;
            msExptAnalyzed(exptN).av(auditoryTrials).cue(valid).highThreshHR = [nan nan];
            msExptAnalyzed(exptN).av(auditoryTrials).cue(invalid).highThreshHR = [nan nan];

        end
        
        if isfield(msCmlvData,'tVisTargets')
            msCmlvData.tVisTargets = cat(2,msCmlvData.tVisTargets,tVisTargets);
            msCmlvData.tAudTargets = cat(2,msCmlvData.tAudTargets,tAudTargets);
            msCmlvData.tInvVisTargets = cat(2,msCmlvData.tInvVisTargets,tInvVisTargets);
            msCmlvData.tInvAudTargets = cat(2,msCmlvData.tInvAudTargets,tInvAudTargets);
            msCmlvData.hit = cat(2,msCmlvData.hit,hit);
            msCmlvData.miss = cat(2,msCmlvData.miss,miss);
            msCmlvData.fa = cat(2,msCmlvData.fa,fa);
            msCmlvData.invHit = cat(2,msCmlvData.invHit,invHit);
            msCmlvData.invMiss = cat(2,msCmlvData.invMiss,invMiss);
            msCmlvData.valRT = cat(2,msCmlvData.valRT,valRT);
            msCmlvData.invRT = cat(2,msCmlvData.invRT,invRT);
            msCmlvData.lastStimRT = cat(2,msCmlvData.lastStimRT,lastStimRT);
            msCmlvData.oneStimBackRT = cat(2,msCmlvData.oneStimBackRT,oneStimBackRT);
            msCmlvData.twoStimBackRT = cat(2,msCmlvData.twoStimBackRT,twoStimBackRT);
         
            msCmlvData.valTargetTimeMs = cat(2,msCmlvData.valTargetTimeMs,targetTimeMs);
            msCmlvData.invTargetTimeMs = cat(2,msCmlvData.invTargetTimeMs,invTargetTimeMs);
            msCmlvData.visMatchedHighThreshID = cat(2,msCmlvData.visMatchedHighThreshID,visMatchedHighThreshID);
            msCmlvData.visMatchedHighThreshOutcome = cat(2,msCmlvData.visMatchedHighThreshOutcome,visMatchedHighThreshOutcome);
            msCmlvData.invVisHighThreshID = cat(2,msCmlvData.invVisHighThreshID,invVisHighThreshID);
            msCmlvData.audMatchedHighThreshID = cat(2,msCmlvData.audMatchedHighThreshID,audMatchedHighThreshID);
            msCmlvData.audMatchedHighThreshOutcome = cat(2,msCmlvData.audMatchedHighThreshOutcome,audMatchedHighThreshOutcome);
            msCmlvData.invAudHighThreshID = cat(2,msCmlvData.invAudHighThreshID,invAudHighThreshID);
        else
            msCmlvData.tVisTargets = tVisTargets;
            msCmlvData.tAudTargets = tAudTargets;
            msCmlvData.tInvVisTargets = tInvVisTargets;
            msCmlvData.tInvAudTargets = tInvAudTargets;
            msCmlvData.hit = hit;
            msCmlvData.miss = miss;
            msCmlvData.fa = fa;
            msCmlvData.invHit = invHit;
            msCmlvData.invMiss = invMiss;
            msCmlvData.valRT = valRT;
            msCmlvData.invRT = invRT;
            msCmlvData.lastStimRT = lastStimRT;
            msCmlvData.oneStimBackRT = oneStimBackRT;
            msCmlvData.twoStimBackRT = twoStimBackRT;
            msCmlvData.valTargetTimeMs = targetTimeMs;
            msCmlvData.invTargetTimeMs = invTargetTimeMs;
            msCmlvData.visMatchedHighThreshID = visMatchedHighThreshID;
            msCmlvData.visMatchedHighThreshOutcome = visMatchedHighThreshOutcome;
            msCmlvData.invVisHighThreshID = invVisHighThreshID;
            msCmlvData.audMatchedHighThreshID = audMatchedHighThreshID;
            msCmlvData.audMatchedHighThreshOutcome = audMatchedHighThreshOutcome;
            msCmlvData.invAudHighThreshID = invAudHighThreshID;
        end
    end
    save(fullfile(fnout,[savedDataName 'Analyzed']),...
        'exptInd', 'msExptAnalyzed', 'msCmlvData')
    %%
    if doPlot
        n = sum(exptInd);
        %%
        setFigParams4Print('portrait')
        set(0,'defaultAxesFontSize',12)
        nb = 5;
        figure
        suptitle(sprintf('MS %s: Thresholds for each expt',ms2analyze{im}))
        dv = [];
        da = [];
        for iexp = 1:n
            dv = cat(2,dv,msExptAnalyzed(iexp).av(visualTrials).threshold);
            da = cat(2,da,msExptAnalyzed(iexp).av(auditoryTrials).threshold);
        end
        subplot 221
        h = histogram(dv,nb);
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        title(sprintf('Visual Trials (%s%s%s)',...
            num2str(round(mean(dv),2,'Significant')),'\pm',num2str(round(std(dv),2,'Significant'))))
        figYAxis([],'N Expt',visThresh_lim)
        figXAxis([],'Expt 50% Threshold',[])
        figAxForm
        hold on
        vline(mean(dv),'k--')
        subplot 222
        h = histogram(da,nb);
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        title(sprintf('Auditory Trials (%s%s%s)',...
            num2str(round(mean(da),2,'Significant')),'\pm',num2str(round(std(da),2,'Significant'))))
        figYAxis([],'N Expt',audThresh_lim)
        figXAxis([],'Expt 50% Threshold',[])
        figAxForm
        hold on
        vline(mean(da),'k--')
        dv = [];
        da = [];
        for iexp = 1:n
            dv = cat(2,dv,msExptAnalyzed(iexp).av(visualTrials).highThreshold);
            da = cat(2,da,msExptAnalyzed(iexp).av(auditoryTrials).highThreshold);
        end
        subplot 223
        h = histogram(dv,nb);
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        title(sprintf('Visual Trials (%s%s%s)',...
            num2str(round(mean(dv),2,'Significant')),'\pm',num2str(round(std(dv),2,'Significant'))))
        figYAxis([],'N Expt',visThresh_lim)
        figXAxis([],sprintf('Expt %s Threshold',num2str(highThreshold*100)),[])
        figAxForm
        hold on
        vline(mean(dv),'k--')
        subplot 224
        h = histogram(da,nb);
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        title(sprintf('Auditory Trials (%s%s%s)',...
            num2str(round(mean(da),2,'Significant')),'\pm',num2str(round(std(da),2,'Significant'))))
        figYAxis([],'N Expt',audThresh_lim)
        figXAxis([],sprintf('Expt %s Threshold',num2str(highThreshold*100)),[])
        figAxForm
        hold on
        vline(mean(da),'k--')

        print(fullfile(fnout,[fileName 'exptThresholds']),'-dpdf','-fillpage')
        %%
        setFigParams4Print('landscape')
        set(0,'defaultAxesFontSize',16)
        figure
        suptitle(sprintf('MS %s: Threshold is %s pct, above (black) or below(blue) for each expt',ms2analyze{im},num2str(highThreshold*100)))
        for i = 1:2
            subplot 121
            x = [];
            y = [];
            for iexp = 1:n
                x = cat(2,x,msExptAnalyzed(iexp).av(visualTrials).cue(valid).highThreshHR(i).*100);
                y = cat(2,y,msExptAnalyzed(iexp).av(visualTrials).cue(invalid).highThreshHR(i).*100);
            end
            h = scatter(x,y,'o');
            hold on
            hh = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'o');
            if i == 1
                h.MarkerEdgeColor = 'b';
                h.MarkerFaceColor = 'b';
                hh.Color = 'b';
                title('Visual Trials')
                figXAxis([],'Valid HR (%)',HR_lim,HR_label)
                figYAxis([],'Invalid HR (%)',HR_lim,HR_label,HR_label)
                figAxForm
            else
                h.MarkerEdgeColor = 'k';
                h.MarkerFaceColor = 'k';
                hh.Color = 'k';
                plot(HR_lim,HR_lim,'k--')
            end
            subplot 122
            x = [];
            y = [];
            for iexp = 1:n
                x = cat(2,x,msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).highThreshHR(i).*100);
                y = cat(2,y,msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).highThreshHR(i).*100);
            end
            h = scatter(x,y,'o');
            hold on
            hh = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'o');
            if i == 1
                h.MarkerEdgeColor = 'b';
                h.MarkerFaceColor = 'b';
                hh.Color = 'b';
                title('Auditory Trials')
                figXAxis([],'Valid HR (%)',HR_lim,HR_label)
                figYAxis([],'Invalid HR (%)',HR_lim,HR_label,HR_label)
                figAxForm
            else
                h.MarkerEdgeColor = 'k';
                h.MarkerFaceColor = 'k';
                hh.Color = 'k';
                plot(HR_lim,HR_lim,'k--')
            end
        end

        print(fullfile(fnout,[fileName 'HReaExpt_highThreshold']),'-dpdf','-fillpage')
        %%    
        setFigParams4Print('landscape')
        set(0,'defaultAxesFontSize',12)
        figure
        xv = [];
        xa = [];
        dv = [];
        da = [];
        for iexp = 1:n
            xv = cat(2,xv,...
                msExptAnalyzed(iexp).av(visualTrials).cue(valid).targets - ...
                msExptAnalyzed(iexp).av(visualTrials).threshold);
            dv = cat(2,dv,msExptAnalyzed(iexp).av(visualTrials).cue(valid).HR.*100);
            xa = cat(2,xa,...
                msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).targets - ...
                msExptAnalyzed(iexp).av(auditoryTrials).threshold);
            da = cat(2,da,msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).HR.*100);
            if ~isempty(msExptAnalyzed(iexp).av(auditoryTrials).threshold)
                xa = cat(2,xa,...
                    msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).targets - ...
                    msExptAnalyzed(iexp).av(auditoryTrials).threshold);
                da = cat(2,da,msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).HR.*100);
            end
        end
        subplot 231
        h = scatter(xv,dv,'o');
        figXAxis([],'Target Diff From Threshold (\circ)',[-90 90])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Visual: Valid')
        hold on
        vline(0,'k--')
        subplot 234
        h = scatter(xa,da,'o');
        figXAxis([],'Target Diff From Threshold (a.u.)',[-0.2 0.2])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Auditory: Valid')
        hold on
        vline(0,'k--')
        xv = [];
        xa = [];
        dv1 = [];
        dv2 = [];
        da1 = [];
        da2 = [];
        for iexp = 1:n
            xv = cat(2,xv,...
                msExptAnalyzed(iexp).av(visualTrials).cue(invalid).targets - ...
                msExptAnalyzed(iexp).av(visualTrials).threshold);
            dv1 = cat(2,dv2,msExptAnalyzed(iexp).av(visualTrials).cue(invalid).matchedHR.*100);
            dv2 = cat(2,dv2,...
                (msExptAnalyzed(iexp).av(visualTrials).cue(valid).matchedHR.*100) - ...
                (msExptAnalyzed(iexp).av(visualTrials).cue(invalid).matchedHR.*100));
            xa = cat(2,xa,...
                msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).targets - ...
                msExptAnalyzed(iexp).av(auditoryTrials).threshold);
            da1 = cat(2,da2,msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).matchedHR.*100);
            da2 = cat(2,da2,...
                (msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).matchedHR.*100) - ...
                (msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).matchedHR.*100));
        end

        subplot 232
        h = scatter(xv,dv1,'o');
        figXAxis([],'Target Diff From Threshold (\circ)',[-90 90])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Visual: Invalid')
        hold on
        vline(0,'k--')
        subplot 233
        h = scatter(xv,dv2,'o');
        figXAxis([],'Target Diff From Threshold (\circ)',[-90 90])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        title('Visual: Valid - Invalid')
        hold on
        vline(0,'k--')
        hline(0,'k--')
        subplot 235
        h = scatter(xa,da1,'o');
        figXAxis([],'Target Diff From Threshold (a.u.)',[-0.2 0.2])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Auditory: Invalid')
        hold on
        vline(0,'k--')
        subplot 236
        h = scatter(xa,da2,'o');
        figXAxis([],'Target Diff From Threshold (a.u.)',[-0.2 0.2])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        title('Auditory: Valid - Invalid')
        hold on
        vline(0,'k--')
        hline(0,'k--')

        print(fullfile(fnout,[fileName 'HRxTargetDiffFromThresh']),'-dpdf','-fillpage')

        %%       
        setFigParams4Print('portrait')
        set(0,'defaultAxesFontSize',12)
        figure
        suptitle('Visual Trials')
        x = [];
        d1 = [];
        d2 = [];
        for iexp = 1:n
            x = cat(2,x,msExptAnalyzed(iexp).av(visualTrials).falseAlarmRate.*100);
            d1 = cat(1,d1,msExptAnalyzed(iexp).av(visualTrials).cue(valid).highThreshHR.*100);
            d2 = cat(1,d2,msExptAnalyzed(iexp).av(visualTrials).cue(invalid).highThreshHR.*100);
        end
        subplot 221
        h = scatter(x,d1(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,d1(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Valid')
        subplot 222
        h = scatter(x,d2(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,d2(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Invalid')
        subplot 223
        y = d1 - d2;    
        h = scatter(x,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')
        subplot 224
        y = x;
        x = [];
        for iexp = 1:n
            x = cat(2,x,...
                (sum(msExptAnalyzed(iexp).av(visualTrials).cue(invalid).nTrials)./...
                (sum(msExptAnalyzed(iexp).av(visualTrials).cue(invalid).nTrials) + ...
                sum(msExptAnalyzed(iexp).av(visualTrials).cue(valid).nTrials))).*100);
        end
        subplot 224   
        h = scatter(x,y,'ko');
        figXAxis([],'% Invalid Trials (of matched trials)',[])
        figYAxis([],'False Alarm Rate (%)',[])
        figAxForm
        title('All Matched Trials')

        print(fullfile(fnout,[fileName 'visFAR']),'-dpdf','-fillpage')

        figure
        suptitle('Auditory Trials')
        x = [];
        d1 = [];
        d2 = [];
        for iexp = 1:n
            x = cat(2,x,msExptAnalyzed(iexp).av(auditoryTrials).falseAlarmRate.*100);
            d1 = cat(1,d1,msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).highThreshHR.*100);
            d2 = cat(1,d2,msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).highThreshHR.*100);
        end
        subplot 221
        h = scatter(x,d1(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,d1(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Valid')
        subplot 222
        h = scatter(x,d2(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,d2(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',HR_lim)
        figAxForm
        title('Invalid')
        subplot 223
        y = d1 - d2;    
        h = scatter(x,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'False Alarm Rate (%)',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')
        subplot 224
        y = x;
        x = [];
        for iexp = 1:n
            x = cat(2,x,...
                (sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).nTrials)./...
                (sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).nTrials) + ...
                sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).nTrials))).*100);
        end
        subplot 224   
        h = scatter(x,y,'ko');
        figXAxis([],'% Invalid Trials (of matched trials)',[])
        figYAxis([],'False Alarm Rate (%)',[])
        figAxForm
        title('All Matched Trials')

        print(fullfile(fnout,[fileName 'audFAR']),'-dpdf','-fillpage')

        %%
        setFigParams4Print('portrait')
        set(0,'defaultAxesFontSize',16)
        figure
        suptitle('Visual Trials')
        x1 = [];
        x2 = [];
        d1 = [];
        d2 = [];
        for iexp = 1:n
            x1 = cat(2,x1,getDateFromNumberForm(msExptAnalyzed(iexp).date));
            x2 = cat(2,x2,...
                (sum(msExptAnalyzed(iexp).av(visualTrials).cue(invalid).nTrials)./...
                (sum(msExptAnalyzed(iexp).av(visualTrials).cue(invalid).nTrials) + ...
                sum(msExptAnalyzed(iexp).av(visualTrials).cue(valid).nTrials))).*100);
            d1 = cat(1,d1,...
                msExptAnalyzed(iexp).av(visualTrials).cue(valid).highThreshHR.*100);
            d2 = cat(1,d2,...
                msExptAnalyzed(iexp).av(visualTrials).cue(invalid).highThreshHR.*100);
        end
        x1 = cumsum([0 days(duration(diff(x1)))]);
        y = d1 - d2;  
        subplot 121
        h = scatter(x1,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x1,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'Days from 1st day of dataset',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')
        subplot 122
        h = scatter(x2,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x2,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'% Invalid Trials in Expt',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')

        print(fullfile(fnout,[fileName 'visHRbyExptDay']),'-dpdf','-fillpage')

        figure
        suptitle('Auditory Trials')
        x1 = [];
        x2 = [];
        d1 = [];
        d2 = [];
        for iexp = 1:n
            x1 = cat(2,x1,getDateFromNumberForm(msExptAnalyzed(iexp).date));
            x2 = cat(2,x2,...
                (sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).nTrials)./...
                (sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).nTrials) + ...
                sum(msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).nTrials))).*100);
            d1 = cat(1,d1,...
                msExptAnalyzed(iexp).av(auditoryTrials).cue(valid).highThreshHR.*100);
            d2 = cat(1,d2,...
                msExptAnalyzed(iexp).av(auditoryTrials).cue(invalid).highThreshHR.*100);
        end
        x1 = cumsum([0 days(duration(diff(x1)))]);
        y = d1 - d2;  
        subplot 121
        h = scatter(x1,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x1,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'Days from 1st day of dataset',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')
        subplot 122
        h = scatter(x2,y(:,1),'o');
        h.MarkerFaceColor = 'b';
        hold on
        h = scatter(x2,y(:,2),'o');
        h.MarkerFaceColor = 'k';
        figXAxis([],'% Invalid Trials in Expt',[])
        figYAxis([],'Hit Rate (%)',subHR_lim)
        figAxForm
        hline(0,'k--')
        title('Valid - Invalid')

        print(fullfile(fnout,[fileName 'audHRbyExptDay']),'-dpdf','-fillpage')

        %%
        visTargets = unique(msCmlvData.tVisTargets);
        visTargets = visTargets(2:end);
        audTargets = unique(msCmlvData.tAudTargets);
        audTargets = audTargets(2:end);

        visLevels_lim = [min(visTargets)-1 110];
        visLevels_label = [10 100];
        audLevels_lim = [min(audTargets(audTargets > 0.00001))-...
            (0.5*min(audTargets(audTargets > 0.00001))) 1.1];
        audLevels_label = [0 0.001 0.01 0.1 1];

        %%
        visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
        audBinEdges = exp(linspace(...
            log(min(audTargets(audTargets > 0.00001))-...
            (0.5*min(audTargets(audTargets > 0.00001)))),...
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
        valVisRT = nan(1,nBins);
        invVisRT = nan(1,nBins);
        valAudRT = nan(1,nBins);
        invAudRT = nan(1,nBins);
        valVisRTSte = nan(1,nBins);
        invVisRTSte = nan(1,nBins);
        valAudRTSte = nan(1,nBins);
        invAudRTSte = nan(1,nBins);

        for ibin = 1:nBins
            ind = (msCmlvData.hit | msCmlvData.miss) & visBinInd == ibin;
            if sum(ind) > minTrN_ms
                nVisHits(ibin) = sum(msCmlvData.hit & visBinInd == ibin);
                nVisMisses(ibin) = sum(msCmlvData.miss & visBinInd == ibin);
                visTargetsBinned(ibin) = mean(msCmlvData.tVisTargets(ind));
                visTargetsSte(ibin) = ste(msCmlvData.tVisTargets(ind),2);
                valVisRT(ibin) = mean(...
                    msCmlvData.valRT(msCmlvData.hit & visBinInd == ibin));
                valVisRTSte(ibin) = ste(...
                    msCmlvData.valRT(msCmlvData.hit & visBinInd == ibin),2);
            end

            ind = (msCmlvData.hit | msCmlvData.miss) & audBinInd == ibin;
            if sum(ind) > minTrN_ms
                nAudHits(ibin) = sum(msCmlvData.hit & audBinInd == ibin);
                nAudMisses(ibin) = sum(msCmlvData.miss & audBinInd == ibin);
                audTargetsBinned(ibin) = mean(msCmlvData.tAudTargets(ind));
                audTargetsSte(ibin) = ste(msCmlvData.tAudTargets(ind),2);
                valAudRT(ibin) = mean(...
                    msCmlvData.valRT(msCmlvData.hit & audBinInd == ibin));
                valAudRTSte(ibin) = ste(...
                    msCmlvData.valRT(msCmlvData.hit & audBinInd == ibin),2);
            end

            ind = (msCmlvData.invHit | msCmlvData.invMiss) & invVisBinInd == ibin;
            if sum(ind) > minTrN_ms
                nInvVisHits(ibin) = sum(invVisBinInd == ibin & msCmlvData.invHit);
                nInvVisMisses(ibin) = sum(invVisBinInd == ibin & msCmlvData.invMiss);
                invVisTargetsBinned(ibin) = mean(msCmlvData.tInvVisTargets(ind));
                invVisTargetsSte(ibin) = ste(msCmlvData.tInvVisTargets(ind),2);
                invVisRT(ibin) = mean(...
                    msCmlvData.invRT(msCmlvData.invHit & invVisBinInd == ibin));
                invVisRTSte(ibin) = ste(...
                    msCmlvData.invRT(msCmlvData.invHit & invVisBinInd == ibin),2);
            end

            ind = (msCmlvData.invHit | msCmlvData.invMiss) & invAudBinInd == ibin;
            if sum(ind) > minTrN_ms  
                nInvAudHits(ibin) = sum(invAudBinInd == ibin & msCmlvData.invHit);
                nInvAudMisses(ibin) = sum(invAudBinInd == ibin & msCmlvData.invMiss);
                invAudTargetsBinned(ibin) = mean(msCmlvData.tInvAudTargets(ind));
                invAudTargetsSte(ibin) = ste(msCmlvData.tInvAudTargets(ind),2);
                invAudRT(ibin) = mean(...
                    msCmlvData.invRT(msCmlvData.invHit & invAudBinInd == ibin));
                invAudRTSte(ibin) = ste(...
                    msCmlvData.invRT(msCmlvData.invHit & invAudBinInd == ibin),2);
            end
        end
        visInd = ~isnan(nVisHits);
        invVisInd = ~isnan(nInvVisHits); 
        audInd = ~isnan(nAudHits);
        invAudInd = ~isnan(nInvAudHits); 

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

        %% 
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

        %%
        visHighThreshDeg = msVisXGrid(find(msVisFit.modelFun(msVisFit.coefEsts,msVisXGrid) ...
            > highThreshold,1));
        audHighThreshAmp = msAudXGrid(find(msAudFit.modelFun(msAudFit.coefEsts,msAudXGrid) ...
            > highThreshold,1));

        %% all trials attention
        [visAttnP, valHR_vis, invHR_vis, valCI_vis, invCI_vis] = ...
            allTrialsAttnTest(msCmlvData.tVisTargets,msCmlvData.tInvVisTargets,...
            msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
        [audAttnP, valHR_aud, invHR_aud, valCI_aud, invCI_aud] = ...
            allTrialsAttnTest(msCmlvData.tAudTargets,msCmlvData.tInvAudTargets,...
            msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);

        figure
        subplot 121
        x = 1:2;
        y = cat(2,valHR_vis,invHR_vis).*100;
        ylerr = y - cat(2,valCI_vis(1),invCI_vis(1)).*100;
        yuerr = cat(2,valCI_vis(2),invCI_vis(2)).*100 - y;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = threshColors{i};
        hold on
        figXAxis([],'',[0 3],x,{'Valid';'Invalid'});
        figYAxis([],'Matched Hit Rate (%)',HR_lim)
        figAxForm([],0)
        title('Visual Trials')
        if visAttnP < 0.05
            text(mean(x), mean(y),sprintf('* p = %s (%s)',...
                num2str(round(visAttnP,3,'significant')),...
                num2str(round(y(1)-y(2),3,'significant'))))
        else
            text(mean(x), mean(y),sprintf('p = %s (%s)',...
                num2str(round(visAttnP,3,'significant')),...
                num2str(round(y(1)-y(2),3,'significant'))))
        end
        disp(sprintf('vis all: p = %s (%s)',...
            num2str(round(visAttnP,3,'significant')),...
            num2str(round(y(1)-y(2),3,'significant'))))
        subplot 122
        x = 1:2;
        y = cat(2,valHR_aud,invHR_aud).*100;
        ylerr = y - cat(2,valCI_aud(1),invCI_aud(1)).*100;
        yuerr = cat(2,valCI_aud(2),invCI_aud(2)).*100 - y;
        h = errorbar(x,y,ylerr,yuerr,'ko-');
        h.MarkerFaceColor = threshColors{i};
        hold on
        figXAxis([],'',[0 3],x,{'Valid';'Invalid'});
        figYAxis([],'Matched Hit Rate (%)',HR_lim)
        figAxForm([],0)
        title('Auditory Trials')
        if audAttnP < 0.05
            text(mean(x), mean(y),sprintf('* p = %s (%s)',...
                num2str(round(audAttnP,3,'significant')),...
                num2str(round(y(1)-y(2),3,'significant'))))
        else
            text(mean(x), mean(y),sprintf('p = %s (%s)',...
                num2str(round(audAttnP,3,'significant')),...
                num2str(round(y(1)-y(2),3,'significant'))))
        end
        disp(sprintf('aud all: p = %s (%s)',...
            num2str(round(audAttnP,3,'significant')),...
            num2str(round(y(1)-y(2),3,'significant'))))
        %%
        visAttnP = belowThreshAttnTest(visHighThreshDeg,...
            msCmlvData.tVisTargets,msCmlvData.tInvVisTargets,...
            msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);
        audAttnP = belowThreshAttnTest(audHighThreshAmp,...
            msCmlvData.tAudTargets,msCmlvData.tInvAudTargets,...
            msCmlvData.hit,msCmlvData.miss,msCmlvData.invHit,msCmlvData.invMiss);


        %%
        nExpAnalyzed = size(msExptAnalyzed,2);
        bxSummaryTable = table(nExpAnalyzed,nansum(nVisHits+nVisMisses),...
            nansum(nInvVisHits+nInvVisMisses),nExpAnalyzed,...
            nansum(nAudHits+nAudMisses),nansum(nInvAudHits+nInvAudMisses));
        bxSummaryTable.Properties.VariableNames = {'nVisSessions','nVisValidIncluded',....
            'nVisInvalidIncl','nAudSessions','nAudValidInlcuded','nAudInvalidIncluded'};

        %%
        setFigParams4Print('landscape')
        set(0,'defaultAxesFontSize',12)
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
        vline(visHighThreshDeg,'k--',sprintf('threshold = %s',num2str(visHighThreshDeg)))
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
        figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
        figAxForm
        title({'Visual Trials';sprintf('attn: p = %s', num2str(visAttnP))})

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
        vline(audHighThreshAmp,'k--',sprintf('threshold = %s',num2str(audHighThreshAmp)))
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
        figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
        figAxForm
        title({'Auditory Trials';sprintf('attn: p = %s', num2str(audAttnP))})

        subplot 133
        ax = gca;
        ax.Visible = 'off';
        h = text(0.6,0.5,bxSummaryTable.Properties.VariableNames',...
            'FontWeight', 'bold', 'FontSize', 10,'HorizontalAlignment','right');
        bxSummaryNs = cellfun(@num2str,...
            {nExpAnalyzed,nansum(nVisHits+nVisMisses),...
            nansum(nInvVisHits+nInvVisMisses),nExpAnalyzed,...
            nansum(nAudHits+nAudMisses),nansum(nInvAudHits+nInvAudMisses)},...
            'unif',0)';
        text(0.7, 0.5,bxSummaryNs,'FontWeight', 'bold', 'FontSize', 10);

        print(fullfile(fnout, [fileName 'allDataHRxDifficulty']),'-dpdf','-fillpage')

        %%
        nMatchedValidHR_vis = nan(1,2);
        nMatchedInvalidHR_vis = nan(1,2);
        nMatchedValidCI_vis = nan(2,2);
        nMatchedInvalidCI_vis = nan(2,2);    
        nMatchedValidHR_aud = nan(1,2);
        nMatchedInvalidHR_aud = nan(1,2);
        nMatchedValidCI_aud = nan(2,2);
        nMatchedInvalidCI_aud = nan(2,2);  
        matchedHighThreshTest_vis = nan(1,2);   
        matchedHighThreshTest_aud = nan(1,2);    
        for i = 1:2
            if i == 1
                visValInd = msCmlvData.visMatchedHighThreshID == 1;
                audValInd = msCmlvData.audMatchedHighThreshID == 1;
                visInvInd = msCmlvData.tInvVisTargets < visHighThreshDeg & ...
                    msCmlvData.tInvVisTargets > 0;
                audInvInd = msCmlvData.tInvAudTargets < audHighThreshAmp & ...
                    msCmlvData.tInvAudTargets > 0;  
            elseif i == 2      
                visValInd = msCmlvData.visMatchedHighThreshID == 2;
                audValInd = msCmlvData.audMatchedHighThreshID == 2;
                visInvInd = msCmlvData.tInvVisTargets > visHighThreshDeg;
                audInvInd = msCmlvData.tInvAudTargets > audHighThreshAmp;  
            end

            nH = sum(msCmlvData.visMatchedHighThreshOutcome(visValInd) == 1);
            nM = sum(msCmlvData.visMatchedHighThreshOutcome(visValInd) == 0);
            [nMatchedValidHR_vis(i), nMatchedValidCI_vis(:,i)] = binofit(nH,nH+nM);

            nH = sum(visInvInd & msCmlvData.invHit);
            nM = sum(visInvInd & msCmlvData.invMiss);
            [nMatchedInvalidHR_vis(i), nMatchedInvalidCI_vis(:,i)] = binofit(nH,nH+nM);

            matchedHighThreshTest_vis(i) = binocdf(nH,nH+nM,nMatchedValidHR_vis(i));

            nH = sum(msCmlvData.audMatchedHighThreshOutcome(audValInd) == 1);
            nM = sum(msCmlvData.audMatchedHighThreshOutcome(audValInd) == 0);
            [nMatchedValidHR_aud(i), nMatchedValidCI_aud(:,i)] = binofit(nH,nH+nM);

            nH = sum(audInvInd & msCmlvData.invHit);
            nM = sum(audInvInd & msCmlvData.invMiss);
            [nMatchedInvalidHR_aud(i), nMatchedInvalidCI_aud(:,i)] = binofit(nH,nH+nM);

            matchedHighThreshTest_aud(i) = binocdf(nH,nH+nM,nMatchedValidHR_aud(i));
        end


        setFigParams4Print('portrait')
        set(0,'defaultAxesFontSize',12)
        figure
        subplot 121
        for i = 1:2
            x = 1:2;
            y = cat(2,nMatchedValidHR_vis(i),nMatchedInvalidHR_vis(i)).*100;
            h = plot(x,y,[threshColors{i} 'o-']);
            h.MarkerFaceColor = threshColors{i};
            hold on
            if matchedHighThreshTest_vis(i) < 0.05
                text(mean(x), mean(y),sprintf('* p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_vis(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            else
                text(mean(x), mean(y),sprintf('p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_vis(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            end
            if i == 1
                disp(sprintf('vis below: p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_vis(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            end
        end
        figXAxis([],'',[0 3],x,{'Valid';'Invalid'});
        figYAxis([],'Matched Hit Rate (%)',HR_lim)
        figAxForm([],0)
        title('Visual Trials')
        subplot 122
        for i = 1:2
            x = 1:2;
            y = cat(2,nMatchedValidHR_aud(i),nMatchedInvalidHR_aud(i)).*100;
            h = plot(x,y,[threshColors{i} 'o-']);
            h.MarkerFaceColor = threshColors{i};
            hold on
            if matchedHighThreshTest_aud(i) < 0.05
                text(mean(x), mean(y),sprintf('* p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_aud(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            else
                text(mean(x), mean(y),sprintf('p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_aud(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            end
            if i == 1
                disp(sprintf('aud below: p = %s (%s)',...
                    num2str(round(matchedHighThreshTest_aud(i),3,'significant')),...
                    num2str(round(y(1)-y(2),3,'significant'))))
            end
        end
        figXAxis([],'',[0 3],x,{'Valid';'Invalid'});
        figYAxis([],'Matched Hit Rate (%)',HR_lim)
        figAxForm([],0)
        title('Auditory Trials')

        print(fullfile(fnout,[fileName 'allDataMatchedHRxCue']),'-dpdf')
        %% reaction times
        
        setFigParams4Print('landscape')
        set(0,'defaultAxesFontSize',12)
        figure
        h = cdfplot(msCmlvData.valRT(msCmlvData.valRT > 0 & msCmlvData.valRT < 500));
        h.Color = cueColors{valid};
        hold on
        h = cdfplot(msCmlvData.invRT(msCmlvData.invRT > 0 & msCmlvData.invRT < 500));
        h.Color = cueColors{invalid};
        h = cdfplot(msCmlvData.lastStimRT(msCmlvData.lastStimRT > 0 & ...
            msCmlvData.lastStimRT < 500 & msCmlvData.fa & ...
             ~(msCmlvData.invHit |msCmlvData.invMiss)));
        h.Color = cueColors{valid};
        h.LineStyle = '--';
        figXAxis([],'Reaction Time (ms)',[0 600])
        figYAxis([],'Fraction of Hits',[0 1])
        figAxForm
        title('All Trials (0 < RT < 500)') 
        legend({'Valid','Invalid'},'location','northeastoutside')        
        print(fullfile(fnout,[fileName 'reactTimeDistributionAllTrials']),'-dpdf')
        
        setFigParams4Print('landscape')
        set(0,'defaultAxesFontSize',12)
        figure
        suptitle({mouseName,' (0 < RT < 540)'})
        subplot 221
        ind = ~msCmlvData.hit & ~msCmlvData.miss & msCmlvData.tAudTargets > 0;
        d = cat(2,msCmlvData.lastStimRT(ind & msCmlvData.lastStimRT > 0),...
            msCmlvData.oneStimBackRT(ind & msCmlvData.oneStimBackRT < 540));
        h = histogram(d,[0:20:540],'Normalization','probability');
        h.FaceColor = [0.75 0.75 0];
        hold on
        ind = msCmlvData.valRT > 0 & msCmlvData.valRT < 540 & msCmlvData.tVisTargets > 0;
        h = histogram(msCmlvData.valRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{valid};
        hold on
        ind = msCmlvData.invRT > 0 & msCmlvData.invRT < 540 & msCmlvData.tInvVisTargets > 0;
        h = histogram(msCmlvData.invRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{invalid};
        figXAxis([],'Reaction Time (ms)',[0 540])
        figYAxis([],'Fraction of Trials',[0 0.2])
        figAxForm([],0)
        title('All Visual Trials (Auditory FA)') 
        legend({'FA','Valid','Invalid'},'location','northwest')  
        
        subplot 222
        ind = ~msCmlvData.hit & ~msCmlvData.miss & msCmlvData.tVisTargets > 0;
        d = cat(2,msCmlvData.lastStimRT(ind & msCmlvData.lastStimRT > 0),...
            msCmlvData.oneStimBackRT(ind & msCmlvData.oneStimBackRT < 540));
        h = histogram(d,[0:20:540],'Normalization','probability');
        h.FaceColor = [0.75 0.75 0];
        hold on
        ind = msCmlvData.valRT > 0 & msCmlvData.valRT < 540 & msCmlvData.tAudTargets > 0;
        h = histogram(msCmlvData.valRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{valid};
        hold on
        ind = msCmlvData.invRT > 0 & msCmlvData.invRT < 540 & msCmlvData.tInvAudTargets > 0;
        h = histogram(msCmlvData.invRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{invalid};
        figXAxis([],'Reaction Time (ms)',[0 540])
        figYAxis([],'Fraction of Trials',[0 0.2])
        figAxForm([],0)
        title('All Auditory Trials (Visual FA)') 
        
        subplot 223
        ind = ~msCmlvData.hit & ~msCmlvData.miss & msCmlvData.tAudTargets > 0;
        d = cat(2,msCmlvData.lastStimRT(ind & msCmlvData.lastStimRT > 0),...
            msCmlvData.oneStimBackRT(ind & msCmlvData.oneStimBackRT < 540));
        h = histogram(d,[0:20:540],'Normalization','probability');
        h.FaceColor = [0.75 0.75 0];
        hold on
        ind = msCmlvData.valRT > 0 & msCmlvData.valRT < 540 & msCmlvData.tVisTargets == 23;
        h = histogram(msCmlvData.valRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{valid};
        hold on
        ind = msCmlvData.invRT > 0 & msCmlvData.invRT < 540 & msCmlvData.tInvVisTargets == 23;
        h = histogram(msCmlvData.invRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{invalid};
        figXAxis([],'Reaction Time (ms)',[0 540])
        figYAxis([],'Fraction of Trials',[0 0.2])
        figAxForm([],0)
        title('23 deg Visual Trials (Auditory FA)')   
        
        subplot 224
        ind = ~msCmlvData.hit & ~msCmlvData.miss & msCmlvData.tAudTargets > 0;
        d = cat(2,msCmlvData.lastStimRT(ind & msCmlvData.lastStimRT > 0),...
            msCmlvData.oneStimBackRT(ind & msCmlvData.oneStimBackRT < 540));
        h = histogram(d,[0:20:540],'Normalization','probability');
        h.FaceColor = [0.75 0.75 0];
        hold on
        ind = msCmlvData.valRT > 0 & msCmlvData.valRT < 540 & msCmlvData.tVisTargets == 90;
        h = histogram(msCmlvData.valRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{valid};
        hold on
        ind = msCmlvData.invRT > 0 & msCmlvData.invRT < 540 & msCmlvData.tInvVisTargets == 90;
        h = histogram(msCmlvData.invRT(ind),[0:20:540],'Normalization','probability');
        h.FaceColor = cueColor{invalid};
        figXAxis([],'Reaction Time (ms)',[0 540])
        figYAxis([],'Fraction of Trials',[0 0.2])
        figAxForm([],0)
        title('90 deg Visual Trials (Auditory FA)')   
        print(fullfile(fnout,[fileName 'reactTimeHistogram']),'-dpdf','-fillpage')


        [valVisRTs_matched, invVisRTs_matched] = getMatchedReactTimes(...
            msCmlvData.tVisTargets,msCmlvData.tInvVisTargets,...
            msCmlvData.valRT,msCmlvData.invRT,...
            msCmlvData.hit,msCmlvData.invHit);
        [valAudRTs_matched, invAudRTs_matched] = getMatchedReactTimes(...
            msCmlvData.tAudTargets,msCmlvData.tInvAudTargets,...
            msCmlvData.valRT,msCmlvData.invRT,...
            msCmlvData.hit,msCmlvData.invHit);
        
        setFigParams4Print('portrait')
        set(0,'defaultAxesFontSize',12)
        figure
        suptitle(mouseName)
        subplot 221
        h = cdfplot(valVisRTs_matched);
        h.Color = cueColors{valid};
        hold on
        h = cdfplot(invVisRTs_matched);
        h.Color = cueColors{invalid};
        figXAxis([],'Reaction Time (ms)',[0 600])
        figYAxis([],'Fraction of Hits',[0 1])
        figAxForm
        title('Visual Trials')        
        subplot 222
        h = cdfplot(valAudRTs_matched);
        h.Color = cueColors{valid};
        hold on
        h = cdfplot(invAudRTs_matched);
        h.Color = cueColors{invalid};
        figXAxis([],'Reaction Time (ms)',[0 600])
        figYAxis([],'Fraction of Hits',[0 1])
        figAxForm
        title('Auditory Trials')
                
        subplot 223
        x = visTargetsBinned(visInd);
        xerr = visTargetsSte(visInd);
        y = valVisRT(visInd);
        yerr = valVisRTSte(visInd);
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'o-');
        h.Color = cueColors{valid};
        h.MarkerFaceColor = [1 1 1];
        hold on
        x = invVisTargetsBinned(invVisInd);
        xerr = invVisTargetsSte(invVisInd);
        y = invVisRT(invVisInd);
        yerr = invVisRTSte(invVisInd);
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'o-');
        h.Color = cueColors{invalid};
        h.MarkerFaceColor = [1 1 1];
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
        figYAxis([],'Reaction Time (ms)',rt_lim)
        figAxForm
        title('Visual Trials')

        subplot 224
        x = audTargetsBinned(audInd);
        xerr = audTargetsSte(audInd);
        y = valAudRT(audInd);
        yerr = valAudRTSte(audInd);
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'o-');
        h.Color = cueColors{valid};
        h.MarkerFaceColor = [1 1 1];
        hold on
        x = invAudTargetsBinned(invAudInd);
        xerr = invAudTargetsSte(invAudInd);
        y = invAudRT(invAudInd);
        yerr = invAudRTSte(invAudInd);
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'o-');
        h.Color = cueColors{invalid};
        h.MarkerFaceColor = [1 1 1];
        ax = gca;
        ax.XScale = 'log';
        figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
        figYAxis([],'Reaction Time (ms)',rt_lim)
        figAxForm
        title('Auditory Trials')

        print(fullfile(fnout,[fileName 'reactTimeByDifficulty']),'-dpdf','-fillpage')

        %%
%         [~,~,valTimeGroup] = histcounts(msCmlvData.valTargetTimeMs,timeBins);
%         [~,~,invTimeGroup] = histcounts(msCmlvData.invTargetTimeMs,timeBins);
%         nTimeBins = length(timeBins)-1;
% 
%         nVis_timeBins = nan(2,nTimeBins);
%         nInvVis_timeBins = nan(2,nTimeBins);
%         nAud_timeBins = nan(2,nTimeBins);
%         nInvAud_timeBins = nan(2,nTimeBins);
%         meanTimePerBin = nan(2,nTimeBins);
%         errTimePerBin = nan(2,nTimeBins);
%         for ibin = 1:nTimeBins
%             invInd = invTimeGroup == ibin & msCmlvData.tInvVisTargets > 0;
%             if sum(invInd) >= minTrN_ms
%                 invTargets = unique(msCmlvData.tInvVisTargets(invInd));
%                 valInd = valTimeGroup == ibin & ismember(msCmlvData.tVisTargets,invTargets);
%                 nVis_timeBins(1,ibin) = sum(msCmlvData.hit(valInd));
%                 nVis_timeBins(2,ibin) = sum(msCmlvData.miss(valInd));
%                 nInvVis_timeBins(1,ibin) = sum(msCmlvData.invHit(invInd));
%                 nInvVis_timeBins(2,ibin) = sum(msCmlvData.invMiss(invInd));
%                 meanTimePerBin(visualTrials,ibin) = nanmean(...
%                     cat(2,msCmlvData.valTargetTimeMs(valInd),...
%                     msCmlvData.invTargetTimeMs(invInd)));
%                 errTimePerBin(visualTrials,ibin) = ste(...
%                     cat(2,msCmlvData.valTargetTimeMs(valInd),...
%                     msCmlvData.invTargetTimeMs(invInd)),2);
%             end
% 
%             invInd = invTimeGroup == ibin & msCmlvData.tInvAudTargets > 0;
%             if sum(invInd) >= minTrN_ms
%                 invTargets = unique(msCmlvData.tInvAudTargets(invInd));
%                 valInd = valTimeGroup == ibin & ismember(msCmlvData.tAudTargets,invTargets);
%                 nAud_timeBins(1,ibin) = sum(msCmlvData.hit(valInd));
%                 nAud_timeBins(2,ibin) = sum(msCmlvData.miss(valInd));
%                 nInvAud_timeBins(1,ibin) = sum(msCmlvData.invHit(invInd));
%                 nInvAud_timeBins(2,ibin) = sum(msCmlvData.invMiss(invInd));
%                 meanTimePerBin(auditoryTrials,ibin) = nanmean(...
%                     cat(2,msCmlvData.valTargetTimeMs(valInd),...
%                     msCmlvData.invTargetTimeMs(invInd)));
%                 errTimePerBin(auditoryTrials,ibin) = ste(...
%                     cat(2,msCmlvData.valTargetTimeMs(valInd),...
%                     msCmlvData.invTargetTimeMs(invInd)),2);
%             end
%         end
% 
%         ind = ~isnan(sum(nVis_timeBins));
%         [visHR_timeBins,visCI_timeBins] = binofit(nVis_timeBins(1,ind),...
%             sum(nVis_timeBins(:,ind),1));
%         [invVisHR_timeBins,invVisCI_timeBins] = binofit(nInvVis_timeBins(1,ind),...
%             sum(nInvVis_timeBins(:,ind),1));
%         if any(ind)
%             [audHR_timeBins,audCI_timeBins] = binofit(nAud_timeBins(1,ind),...
%                 sum(nAud_timeBins(:,ind),1));
%             [invAudHR_timeBins,invAudCI_timeBins] = binofit(nInvAud_timeBins(1,ind),...
%                 sum(nInvAud_timeBins(:,ind),1));
%         else
%             audHR_timeBins = nan(1,2);
%             audCI_timeBins = nan(2,2);
%             invAudHR_timeBins = nan(1,2);
%             invAudCI_timeBins = nan(2,2);
%         end
% 
%         setFigParams4Print('landscape')
%         set(0,'defaultAxesFontSize',12)
%         figure
%         suptitle(mouseName)
%         subplot 121
%         ind = ~isnan(meanTimePerBin(visualTrials,:));
%         x = meanTimePerBin(visualTrials,ind);
%         xerr = errTimePerBin(visualTrials,ind);
%         y = visHR_timeBins.*100;
%         ylerr = y - (visCI_timeBins(:,1)'.*100);
%         yuerr = (visCI_timeBins(:,2)'.*100) - y;
%         h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'o-');
%         h.Color = cueColors{valid};
%         h.MarkerFaceColor = [1 1 1];
%         hold on
%         y = invVisHR_timeBins.*100;
%         ylerr = y - (invVisCI_timeBins(:,1)'.*100);
%         yuerr = (invVisCI_timeBins(:,2)'.*100) - y;
%         h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'o-');
%         h.Color = cueColors{invalid};
%         h.MarkerFaceColor = [1 1 1];
%         figXAxis([],'Target Time From Start (ms)',timeLim,timeBins,timeBins)
%         figYAxis([],'HR (%)',HR_lim)
%         figAxForm
%         title('Visual Trials')
%         subplot 122
%         ind = ~isnan(meanTimePerBin(auditoryTrials,:));
%         x = meanTimePerBin(auditoryTrials,ind);
%         xerr = errTimePerBin(auditoryTrials,ind);
%         y = audHR_timeBins.*100;
%         ylerr = y - (audCI_timeBins(:,1)'.*100);
%         yuerr = (audCI_timeBins(:,2)'.*100) - y;
%         h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'o-');
%         h.Color = cueColors{valid};
%         h.MarkerFaceColor = [1 1 1];
%         hold on
%         y = invAudHR_timeBins.*100;
%         ylerr = y - (invAudCI_timeBins(:,1)'.*100);
%         yuerr = (invAudCI_timeBins(:,2)'.*100) - y;
%         h = errorbar(x,y,ylerr,yuerr,xerr,xerr,'o-');
%         h.Color = cueColors{invalid};
%         h.MarkerFaceColor = [1 1 1];
%         figXAxis([],'Target Time From Start (ms)',timeLim,timeBins,timeBins)
%         figYAxis([],'HR (%)',HR_lim)
%         figAxForm
%         title('Auditory Trials')
% 
%         print(fullfile(fnout,[fileName 'HRbyTrialLength']),'-dpdf')
    end
end