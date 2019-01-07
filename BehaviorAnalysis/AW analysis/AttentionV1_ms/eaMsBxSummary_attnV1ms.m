clear all
close all

doLoadPreviousDataset = false;
doCheck4NewDates = true;
doRewarded = true;


bxParams_FSAV_attnV1ms
rc = behavConstsAV;
exptSummaryDir = fullfile(rc.ashley,...
    'Manuscripts','Attention V1','Mouse Info.xlsx');
exptSummaryInfo = readtable(exptSummaryDir);
ms2analyze = exptSummaryInfo.SubjectNumber';


%%
nmice = length(ms2analyze);

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
                save(fullfile(fnout,savedDataName),'msExptInfo')
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
        % get visual trial lapse rates
        if length(msExptInfo(iexp).visHR) > 1
            visHRCutoffPass = msExptInfo(iexp).visHR(end-1:end) < lapse_cutoff;
        else
            visHRCutoffPass = msExptInfo(iexp).visHR < lapse_cutoff;
        end
        % get auditory trial lapse rates
        if length(msExptInfo(iexp).audHR) > 1
            audHRCutoffPass = msExptInfo(iexp).audHR(end-1:end) < lapse_cutoff;
        else
            audHRCutoffPass = msExptInfo(iexp).audHR < lapse_cutoff;
        end
        % is it a catch trial experiment? and does it pass the lapse and
        % early rates?
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
        msExptAnalyzed(exptN).av(visualTrials).nFA = [...
            sum(fa(trueFAInd & visTrials & ~invTrials)),...
            nVisFATrials];
        msExptAnalyzed(exptN).av(auditoryTrials).falseAlarmRate = ...
            sum(fa(trueFAInd & audTrials & ~invTrials))/...
            sum(nFACyc(longEnoughTrials & audTrials & ~invTrials));
        nAudFATrials = sum(nFACyc(longEnoughTrials & audTrials & ~invTrials));
        msExptAnalyzed(exptN).av(auditoryTrials).nFA = [...
            sum(fa(trueFAInd & audTrials & ~invTrials)),...
            nAudFATrials];
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
        
        [visHR_highThreshold, invVisHR_highThreshold,~,~,~,~,...
            highThreshExpt,visMatchHighThreshTrialInd] = ...
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
            [audHR_highThreshold, invAudHR_highThreshold, ~,~,~,~,...
                highThreshExpt ,audMatchHighThreshTrialInd] = ...
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
    
    
    nFA_vis = zeros(1,2);
    nFA_aud = zeros(1,2);
    for iexp = 1:size(msExptAnalyzed,2)
        nFA_vis = nFA_vis + msExptAnalyzed(iexp).av(visualTrials).nFA;
        nFA_aud = nFA_aud + msExptAnalyzed(iexp).av(auditoryTrials).nFA;
    end
    msCmlvData.visNFAandDistractors = nFA_vis;
    msCmlvData.audNFAandDistractors = nFA_aud;
    save(fullfile(fnout,[savedDataName 'Analyzed_attnV1ms']),...
        'exptInd', 'msExptAnalyzed', 'msCmlvData')
end