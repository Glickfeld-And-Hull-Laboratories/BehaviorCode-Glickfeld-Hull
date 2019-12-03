clear all
close all
useRandSeed = false;
ds = 'FSAV_attentionV1';
eval(ds)
rc = behavConstsAV;
imgParams_FSAV

titleStr = ds(6:end);
fnin = fullfile(rc.caOutputDir, ds,[titleStr '_eye_']);
load([fnin 'eyeStruct'])

fnout = fullfile(rc.ashley,'Manuscripts','Attention V1','Matlab Figs');
if useRandSeed
    load(fullfile(fnout,'eyeStats'))
    rng(eyeStats.randGeneratorSeed);
end
%%
mice = unique({msEye.name});
nMice = length(mice);

preEventFr = round(preEventMs_eye*frameRateHz/1000);

%% align experiment data from each mouse
posNames = {'radius','horiz pos','vert pos'};

tcShortStart = cell(3,2,nMice);
tcLongStart = cell(3,2,nMice);
tcTarget = cell(3,2,nMice);

tcShortStart_allTrials = cell(3,nMice);
tcLongStart_allTrials = cell(3,nMice);
tcTarget_allTrials = cell(3,nMice);

for im = 1:nMice
    d = msEye(im);
    fprintf('aligning mouse %s\n',d.name);
    nexp = size(d.expt,2);
    
    start_short_temp = cell(3,2);
    start_long_temp = cell(3,2);
    target_temp = cell(3,2);
    
    start_short_alltemp = cell(3,1);
    start_long_alltemp = cell(3,1);
    target_alltemp = cell(3,1);
    
    for iexp = 1:nexp
        de = d.expt(iexp);
        for ipos = 1:3
            for iav = 1:2
                dp = de.pos(ipos).av(iav);
                start_short_temp{ipos,iav} = cat(2,start_short_temp{ipos,iav},...
                    dp.align(1).shortTC);
                start_long_temp{ipos,iav} = cat(2,start_long_temp{ipos,iav},...
                    dp.align(1).longTC);
                target_temp{ipos,iav} = cat(2,target_temp{ipos,iav},...
                    dp.align(2).TC);
                
                start_short_alltemp{ipos} = cat(2,start_short_alltemp{ipos},...
                    dp.align(1).shortTC);
                start_long_alltemp{ipos} = cat(2,start_long_alltemp{ipos},...
                    dp.align(1).longTC);
                target_alltemp{ipos} = cat(2,target_alltemp{ipos},...
                    dp.align(2).TC);
            end
        end
    end
    tcShortStart(:,:,im) = start_short_temp;
    tcLongStart(:,:,im) = start_long_temp;
    tcTarget(:,:,im) = target_temp;
    tcShortStart_allTrials(:,im) = start_short_alltemp;
    tcLongStart_allTrials(:,im) = start_long_alltemp;
    tcTarget_allTrials(:,im) = target_alltemp;
end

%% timing and plotting params

nFr_short = size(tcShortStart{1,1,1},1);
nFr_long = size(tcLongStart{1,1,1},1);

ttMs_short = ((1:nFr_short) - (preEventFr+nVisDelayFr))./frameRateHz.*1000;
ttMs_long = ((1:nFr_long) - (preEventFr+nVisDelayFr))./frameRateHz.*1000;

posTCLim = [-0.07 0.07];
sizeTCLim = [0.9 1.1];
posLim = [-0.05 0.05];
sizeLim = [0.94 1.06];
%% plot timecourses for each mouse
setFigParams4Print('landscape')
tcAllMice_short = cell(3,2);
tcAllMice_short(:) = deal({nan(length(ttMs_short),nMice)});
tcAllMice_long = cell(3,2);
tcAllMice_long(:) = deal({nan(length(ttMs_long),nMice)});
tcAllMice_target = cell(3,2);
tcAllMice_target(:) = deal({nan(length(ttMs_short),nMice)});
tcAllMice_short_allTrials = cell(3,1);
tcAllMice_short_allTrials(:) = deal({nan(length(ttMs_short),nMice)});
tcAllMice_long_allTrials = cell(3,1);
tcAllMice_long_allTrials(:) = deal({nan(length(ttMs_long),nMice)});
tcAllMice_target_allTrials = cell(3,1);
tcAllMice_target_allTrials(:) = deal({nan(length(ttMs_short),nMice)});
for im = 1:nMice
    eyeFig_AV = figure;
    eyeFig_all = figure;
    for ipos = 1:3
        figure(eyeFig_AV)
        subplot(3,3,ipos)
        x = ttMs_short;
        for iav = 1:2
            y = nanmean(tcShortStart{ipos,iav,im},2);
            yerr = ste(tcShortStart{ipos,iav,im},2);
            hold on
            h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{iav});
            tcAllMice_short{ipos,iav}(:,im) = y;
        end
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Start Align, Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeRespWinFr(1)) x(eyeRespWinFr(end))],'k:')
        figAxForm
        
        figure(eyeFig_all)
        subplot(3,3,ipos)
        x = ttMs_short;
        y = nanmean(tcShortStart_allTrials{ipos,im},2);
        yerr = ste(tcShortStart_allTrials{ipos,im},2);
        hold on
        h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{1});
        tcAllMice_short_allTrials{ipos}(:,im) = y;
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Start Align, Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeRespWinFr(1)) x(eyeRespWinFr(end))],'k:')
        figAxForm
    end

    for ipos = 1:3
        figure(eyeFig_AV)
        subplot(3,3,ipos+3)
        x = ttMs_long;
        for iav = 1:2
            y = nanmean(tcLongStart{ipos,iav,im},2);
            yerr = ste(tcLongStart{ipos,iav,im},2);
            hold on
            h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{iav});
            tcAllMice_long{ipos,iav}(:,im) = y;
        end
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Start Align (long trials), Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeLateRespWinFr(1)) x(eyeLateRespWinFr(end))],'k:')
        figAxForm
        
        figure(eyeFig_all)
        subplot(3,3,ipos+3)
        x = ttMs_long;
        for iav = 1:2
            y = nanmean(tcLongStart_allTrials{ipos,im},2);
            yerr = ste(tcLongStart_allTrials{ipos,im},2);
            hold on
            h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{1});
            tcAllMice_long_allTrials{ipos}(:,im) = y;
        end
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Start Align (long trials), Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeLateRespWinFr(1)) x(eyeLateRespWinFr(end))],'k:')
        figAxForm
    end

    for ipos = 1:3
        figure(eyeFig_AV)
        subplot(3,3,ipos+6)
        x = ttMs_short;
        for iav = 1:2
            y = nanmean(tcTarget{ipos,iav,im},2);
            yerr = ste(tcTarget{ipos,iav,im},2);
            hold on
            h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{iav});
            tcAllMice_target{ipos,iav}(:,im) = y;
        end
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Target Align, Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeRespWinFr(1)) x(eyeRespWinFr(end))],'k:')
        figAxForm
        
        figure(eyeFig_all)
        subplot(3,3,ipos+6)
        x = ttMs_short;
        for iav = 1:2
            y = nanmean(tcTarget_allTrials{ipos,im},2);
            yerr = ste(tcTarget_allTrials{ipos,im},2);
            hold on
            h = shadedErrorBar_chooseColor(x,y,yerr,cueColor{1});
            tcAllMice_target_allTrials{ipos}(:,im) = y;
        end
        figXAxis([],'Time from Start (ms)',[x(1) x(end)])
        if ipos == 1
            figYAxis([],'Norm. Size (mm)',sizeTCLim)
            title({sprintf('Target Align, Mouse %s',mice{im});...
                posNames{ipos}})
        else
            figYAxis([],'\Delta Pos. (mm)',posTCLim)
            title(posNames{ipos})
        end
        vline([x(eyeRespWinFr(1)) x(eyeRespWinFr(end))],'k:')
        figAxForm
    end
    figure(eyeFig_AV)
    print(fullfile(fnout,...
        ['eyeSummary_exampleMouse_AV_' mice{im}]),'-dpdf','-fillpage')
    figure(eyeFig_all)
    print(fullfile(fnout,...
        ['eyeSummary_exampleMouse_allTrials_' mice{im}]),'-dpdf','-fillpage')
end
%% quantify positions for each mouse
earlyAntiResp = cellfun(@(x) mean(x(eyeRespWinFr,:),1),...
    tcAllMice_short,'unif',0);
earlyAntiRespErr = cellfun(@(x) ste(x(eyeRespWinFr,:),1),...
    tcAllMice_short,'unif',0);

lateAntiResp = cellfun(@(x) mean(x(eyeLateRespWinFr,:),1),...
    tcAllMice_long,'unif',0);
lateAntiRespErr = cellfun(@(x) ste(x(eyeLateRespWinFr,:),1),...
    tcAllMice_long,'unif',0);

targetResp = cellfun(@(x) mean(x(eyeRespWinFr,:),1),...
    tcAllMice_target,'unif',0);
targetRespErr = cellfun(@(x) ste(x(eyeRespWinFr,:),1),...
    tcAllMice_target,'unif',0);

earlyAntiResp_allTrials = cellfun(@(x) mean(x(eyeRespWinFr,:),1),...
    tcAllMice_short_allTrials,'unif',0);
earlyAntiRespBL_allTrials = cellfun(@(x) mean(x(eyeBLFr,:),1),...
    tcAllMice_short_allTrials,'unif',0);
[~,earlyTest] = cellfun(@(x,y) ttest(x,y),earlyAntiResp_allTrials,...
    earlyAntiRespBL_allTrials,'unif',0);

lateAntiResp_allTrials = cellfun(@(x) mean(x(eyeLateRespWinFr,:),1),...
    tcAllMice_long_allTrials,'unif',0);
lateAntiRespBL_allTrials = cellfun(@(x) mean(x(eyeBLFr,:),1),...
    tcAllMice_long_allTrials,'unif',0);
[~,lateTest] = cellfun(@(x,y) ttest(x,y),lateAntiResp_allTrials,...
    lateAntiRespBL_allTrials,'unif',0);

targetResp_allTrials = cellfun(@(x) mean(x(eyeRespWinFr,:),1),...
    tcAllMice_target_allTrials,'unif',0);
targetRespBL_allTrials = cellfun(@(x) mean(x(eyeBLFr,:),1),...
    tcAllMice_target_allTrials,'unif',0);
[~,targetTest] = cellfun(@(x,y) ttest(x,y),targetResp_allTrials,...
    targetRespBL_allTrials,'unif',0);

pupilTestColNames = {'Early','Late','Target'};

pupilTest = cell2mat(cat(2,earlyTest,lateTest,targetTest));
%% plot responses
setFigParams4Print('landscape')
figure
suptitle('All Mice, Summary of Eye Responses')
for ipos = 1:3
    subplot (3,3,ipos)
    for iav = 1:2
        x = earlyAntiResp{ipos,visualTrials};
        xerr = earlyAntiRespErr{ipos,visualTrials};
        y = earlyAntiResp{ipos,auditoryTrials};
        yerr = earlyAntiRespErr{ipos,auditoryTrials};
        hold on
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'.');
        xall = mean(x);
        xallerr = ste(x,2);
        yall = mean(y);
        yallerr = ste(y,2);
        h = errorbar(xall,yall,yallerr,yallerr,xallerr,xallerr,'.');
    end
    if ipos == 1        
        figXAxis([],'Visual Trials',sizeLim)
        figYAxis([],'Auditory Trials',sizeLim)
        plot(sizeLim,sizeLim,'k--')
        title([posNames{ipos} ' (Norm)'])
        title({sprintf('Start Align Eye Responses (%s-%s ms)',...
            num2str(round(ttMs_short(eyeRespWinFr(1)))),...
            num2str(round(ttMs_short(eyeRespWinFr(end)))));...
            [posNames{ipos} ' (Norm)']})
    else
        figXAxis([],'Visual Trials',posLim)
        figYAxis([],'Auditory Trials',posLim)
        plot(posLim,posLim,'k--')
        title([posNames{ipos} ' (\Delta in mm)'])
    end
    figAxForm
end

for ipos = 1:3
    subplot (3,3,ipos+3)
    for iav = 1:2
        x = lateAntiResp{ipos,visualTrials};
        xerr = lateAntiRespErr{ipos,visualTrials};
        y = lateAntiResp{ipos,auditoryTrials};
        yerr = lateAntiRespErr{ipos,auditoryTrials};
        hold on
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'.');
        xall = mean(x);
        xallerr = ste(x,2);
        yall = mean(y);
        yallerr = ste(y,2);
        h = errorbar(xall,yall,yallerr,yallerr,xallerr,xallerr,'.');
    end
    if ipos == 1        
        figXAxis([],'Visual Trials',sizeLim)
        figYAxis([],'Auditory Trials',sizeLim)
        plot(sizeLim,sizeLim,'k--')
        title({sprintf('Start Align Eye Responses (%s-%s ms)',...
            num2str(round(ttMs_long(eyeLateRespWinFr(1)))),...
            num2str(round(ttMs_long(eyeLateRespWinFr(end)))));...
            [posNames{ipos} ' (Norm)']})
    else
        figXAxis([],'Visual Trials',posLim)
        figYAxis([],'Auditory Trials',posLim)
        plot(posLim,posLim,'k--')
        title([posNames{ipos} ' (\Delta in mm)'])
    end
    figAxForm
end

for ipos = 1:3
    subplot (3,3,ipos+6)
    for iav = 1:2
        x = targetResp{ipos,visualTrials};
        xerr = targetRespErr{ipos,visualTrials};
        y = targetResp{ipos,auditoryTrials};
        yerr = targetRespErr{ipos,auditoryTrials};
        hold on
        h = errorbar(x,y,yerr,yerr,xerr,xerr,'.');
        xall = mean(x);
        xallerr = ste(x,2);
        yall = mean(y);
        yallerr = ste(y,2);
        h = errorbar(xall,yall,yallerr,yallerr,xallerr,xallerr,'.');
    end
    if ipos == 1        
        figXAxis([],'Visual Trials',sizeLim)
        figYAxis([],'Auditory Trials',sizeLim)
        plot(sizeLim,sizeLim,'k--')
        title({sprintf('Target Align Eye Responses (%s-%s ms)',...
            num2str(round(ttMs_short(eyeRespWinFr(1)))),...
            num2str(round(ttMs_short(eyeRespWinFr(end)))));...
            [posNames{ipos} ' (Norm)']})
    else
%         figXAxis([],'Visual Trials',posLim)
%         figYAxis([],'Auditory Trials',posLim)
        plot(posLim,posLim,'k--')
        title([posNames{ipos} ' (\Delta in mm)'])
    end
    figAxForm
end

print(fullfile(fnout,'eyeSummary_AV_allAttnMice'),'-dpdf','-fillpage')

figure
for ipos = 1:3
    if ipos == 1
        ytest = sizeLim(end);
    else
        ytest = posLim(end);
    end
    subplot(1,3,ipos)
    x = ones(1,nMice);
    y = earlyAntiResp_allTrials{ipos};
    hold on
    plot(x,y,'k.')
    errorbar(1,mean(y),ste(y,2),'.')
    if pupilTest(ipos,1) < eyeAlpha
        plot(1,ytest,'k*')
    end
    y = lateAntiResp_allTrials{ipos};
    plot(x+1,y,'k.')
    errorbar(2,mean(y),ste(y,2),'.')
    if pupilTest(ipos,2) < eyeAlpha
        plot(2,ytest,'k*')
    end
    y = targetResp_allTrials{ipos};
    plot(x+2,y,'k.')
    errorbar(3,mean(y),ste(y,2),'.')
    if pupilTest(ipos,3) < eyeAlpha
        plot(3,ytest,'k*')
    end
    title(posNames{ipos})
    figXAxis([],'',[0 4],1:3,{'Early';'Late';'Target'})
    if ipos == 1
        figYAxis([],'Norm Resp',sizeLim)
        figAxForm
        hline(1,'k--')
    else
        figYAxis([],'\Delta Pos (mm)',posLim)
        figAxForm
        hline(0,'k--')
    end
end
%% save eye stats
eyeStats = struct;
eyeStats.randGeneratorSeed = rng;

for ipos = 1:3
    eyeStats.pos(ipos).name = posNames{ipos};
    eyeStats.pos(ipos).earlyAntiRespDiff = ...
        earlyAntiResp{ipos,visualTrials} - earlyAntiResp{ipos,auditoryTrials};
    [~,eyeStats.pos(ipos).earlyAntiRespTest] = ttest(...
        eyeStats.pos(ipos).earlyAntiRespDiff);
    eyeStats.pos(ipos).lateAntiRespDiff = ...
        lateAntiResp{ipos,visualTrials} - lateAntiResp{ipos,auditoryTrials};
    [~,eyeStats.pos(ipos).lateAntiRespTest] = ttest(...
        eyeStats.pos(ipos).lateAntiRespDiff);
    eyeStats.pos(ipos).targetRespDiff = ...
        targetResp{ipos,visualTrials} - targetResp{ipos,auditoryTrials};
    [~,eyeStats.pos(ipos).targetRespTest] = ttest(...
        eyeStats.pos(ipos).targetRespDiff);
end
visDist = sqrt(sum(cat(1,lateAntiResp{2,visualTrials}.^2,lateAntiResp{3,visualTrials}.^2),1));
audDist = sqrt(sum(cat(1,lateAntiResp{2,auditoryTrials}.^2,lateAntiResp{3,auditoryTrials}.^2),1));
[~,eyeStats.distanceLateAntiTest] = ttest(visDist,audDist);
fprintf('Late Anti: Vis Dist=%s; Aud Dist=%s;p=%s\n',num2str(mean(visDist)),...
    num2str(mean(audDist)),num2str(eyeStats.distanceLateAntiTest))
fprintf('LateAnti Dist Range (um): %s-%s\n',num2str(round(min(cat(2,visDist,audDist)),2,'significant')),...
    num2str(round(max(cat(2,visDist,audDist)),2,'significant')))
visDist = sqrt(sum(cat(1,targetResp{2,visualTrials}.^2,targetResp{3,visualTrials}.^2),1));
audDist = sqrt(sum(cat(1,targetResp{2,auditoryTrials}.^2,targetResp{3,auditoryTrials}.^2),1));
[~,eyeStats.distanceTargetTest] = ttest(visDist,audDist);
fprintf('Target: Vis Dist=%s; Aud Dist=%s;p=%s\n',num2str(mean(visDist)),...
    num2str(mean(audDist)),num2str(eyeStats.distanceTargetTest))
fprintf('Target Dist Range (um): %s-%s\n',num2str(round(min(cat(2,visDist,audDist)),2,'significant')),...
    num2str(round(max(cat(2,visDist,audDist)),2,'significant')))

maxPosChange= nan(1,3);
maxPosChange(2) = max(cat(2,lateAntiResp{2,visualTrials},lateAntiResp{2,auditoryTrials},...
    targetResp{2,visualTrials},targetResp{2,auditoryTrials}));
maxPosChange(3) = max(cat(2,lateAntiResp{3,visualTrials},lateAntiResp{3,auditoryTrials},...
    targetResp{3,visualTrials},targetResp{3,auditoryTrials},...
    targetResp{3,visualTrials},targetResp{3,auditoryTrials}));

for ipos = 1:3
    fprintf('Mean/Err diff early anti %s: %s/%s; p=%s\n',eyeStats.pos(ipos).name,...
        num2str(round(mean(eyeStats.pos(ipos).earlyAntiRespDiff,2),2,'significant')),...
        num2str(round(ste(eyeStats.pos(ipos).earlyAntiRespDiff,2),2,'significant')),...
        num2str(round(eyeStats.pos(ipos).earlyAntiRespTest,2,'significant')))
    fprintf('Mean/Err diff late anti %s: %s/%s; p=%s\n',eyeStats.pos(ipos).name,...
        num2str(round(mean(eyeStats.pos(ipos).lateAntiRespDiff,2),2,'significant')),...
        num2str(round(ste(eyeStats.pos(ipos).lateAntiRespDiff,2),2,'significant')),...
        num2str(round(eyeStats.pos(ipos).lateAntiRespTest,2,'significant')))
    fprintf('Mean/Err diff target %s: %s/%s; p=%s\n',eyeStats.pos(ipos).name,...
        num2str(round(mean(eyeStats.pos(ipos).targetRespDiff,2),2,'significant')),...
        num2str(round(ste(eyeStats.pos(ipos).targetRespDiff,2),2,'significant')),...
        num2str(round(eyeStats.pos(ipos).targetRespTest,2,'significant')))
end

save(fullfile(fnout,'eyeStats'),'eyeStats')