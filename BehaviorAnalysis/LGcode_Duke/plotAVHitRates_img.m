clear all
close all
earlyLT50 = 0;
lapseLT10 = 0;
amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename_fsavImg, [], rc.indexTextCols);
av = behavParamsAV;
mice = unique(xd.Subject);
nMice = length(mice);

% fn = fullfile(rc.fitOutputSummary, [date '_CatchSummary.mat']);
fn = fullfile(rc.fitOutputSummary, ['27-Feb-2017_CatchSummary_fsavImg.mat']);
load(fn)

AFig = figure;
VFig = figure;

colors_mice = parula(nMice);
mice_hr_fig = figure;
easy_all = [];
hard_all = [];
expt_names_all = [];

for imouse = 1:nMice;
    mouse_name = mice(imouse);
    ms_ind = find(mice == mice(imouse));
    if earlyLT50 
        early_ind = find(mouse(imouse).early_mat<0.5);
    else
        early_ind = 1:size(mouse(imouse).early_mat,1);
    end
    if lapseLT10 
        lapse_ind = find(mouse(imouse).HR_ori_mat>0.9);
    else
        lapse_ind = 1:size(mouse(imouse).HR_ori_mat,1);
    end
    
    use_ind = intersect(early_ind,lapse_ind);
    if isempty(use_ind)
        continue 
    end
    disp(size(use_ind))
    input = mouse(imouse).input(use_ind);
    input = concatenateStructures(input');
    
    missedIx = strcmp(input.trialOutcomeCell, 'ignore');
    successIx = strcmp(input.trialOutcomeCell, 'success');
    earliesIx = strcmp(input.trialOutcomeCell, 'failure');
    FAIx = strcmp(input.catchTrialOutcomeCell, 'FA');
    CRIx = strcmp(input.catchTrialOutcomeCell, 'CR');
    earliesIx(FAIx) = 0;
    earliesIx(CRIx) = 0;

    catchDirectionDeg = chop(double(celleqel2mat_padded(input.tCatchGratingDirectionDeg)),2);
%     catchAmplitude = chop(double(celleqel2mat_padded(input.tSoundCatchAmplitude, NaN, 'double')),2);
    catchOris = unique(catchDirectionDeg(~isnan(catchDirectionDeg)));
%     catchAmps = unique(catchAmplitude(~isnan(catchAmplitude)));
    targetDirectionDeg = chop(double(celleqel2mat_padded(input.tGratingDirectionDeg)),2);
%     targetAmplitude = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude, NaN, 'double')),2);
    targetOris = unique(targetDirectionDeg);
%     targetAmps = unique(targetAmplitude(~isnan(targetAmplitude)));
    
% [easy_hr_expt, hard_hr_expt, expt_names] = plotAVHR_expt(mouse,imouse,use_ind,mice_hr_fig,colors_mice,mouse_name);
% 
% easy_all = cat(2,easy_all,easy_hr_expt);
% hard_all = cat(2,hard_all,hard_hr_expt);
% expt_names_all = cat(2,expt_names_all,expt_names);



    cycTimeMs = unique(double(input.stimOffTimeMs)+double(input.stimOnTimeMs));
    nCyc{imouse} = cell2mat(input.tCyclesOn);
    if length(cycTimeMs) > 1
        trialTimeMs{imouse} = bsxfun(@times,double(nCyc{imouse}),double(input.cycTimeMsAllTr));%NaN(size(nCyc{imouse}));
        cycTimeMs_ms{imouse} = double(input.cycTimeMsAllTr);
    else
        trialTimeMs{imouse} = double(cycTimeMs*nCyc{imouse});
        lastStimTimeMs{imouse} = double(cycTimeMs*(nCyc{imouse}-1));
        cycTimeMs_ms{imouse} = repmat(cycTimeMs,1,length(successIx));
    end
      
    tarDeg{imouse} = targetDirectionDeg;
%     tarAmp{imouse} = targetAmplitude;
    cDeg{imouse} = catchDirectionDeg;
%     cAmp{imouse} = catchAmplitude;
    
    sIx_all{imouse} = successIx;
    mIx_all{imouse} = missedIx;
    eIx_all{imouse} = earliesIx;
    faIx_all{imouse} = FAIx;
    crIx_all{imouse} = CRIx;


%     [h_amp, bin_amp] = histc(targetAmplitude, amp_edges);
%     [h_ampc, bin_ampc] = histc(catchAmplitude, amp_edges);
    [h_ori, bin_ori] = histc(targetDirectionDeg, ori_edges);
    [h_oric, bin_oric] = histc(catchDirectionDeg, ori_edges);

    ori_bins = unique(bin_ori);
    avg_ori = zeros(1,max(ori_bins,[],2));
    sem_ori = zeros(1,max(ori_bins,[],2));
    Hits_ori = zeros(1,max(ori_bins,[],2));
    Misses_ori = zeros(1,max(ori_bins,[],2));
    avg_oric = zeros(1,max(ori_bins,[],2));
    sem_oric = zeros(1,max(ori_bins,[],2));
    FAs_ori = zeros(1,max(ori_bins,[],2));
    CRs_ori = zeros(1,max(ori_bins,[],2));
    rct_ori = zeros(1,max(ori_bins,[],2));
    rctsem_ori = zeros(1,max(ori_bins,[],2));
    rct_oric = zeros(1,max(ori_bins,[],2));
    rctsem_oric = zeros(1,max(ori_bins,[],2));
    n_ori = zeros(2,max(ori_bins,[],2));
    for ibin = 1:max(ori_bins,[],2)
        ind = find(bin_ori == ibin);
        if length(ind)>=1
            avg_ori(ibin) = mean(targetDirectionDeg(:,ind),2);
            sem_ori(ibin) = std(targetDirectionDeg,[],2)./sqrt(length(ind));
            Hits_ori(ibin) = sum(successIx(ind),2);
            Misses_ori(ibin) = sum(missedIx(ind),2);
%             rct_ori(ibin) = mean(targetReact(successIx(ind)),2)
%             rctsem_ori(ibin) = std(targetReact(successIx(ind)),[],2)./length(ind);
        end
        indc = find(bin_oric == ibin);
        if length(indc)>=1
            avg_oric(ibin) = mean(catchDirectionDeg(:,indc),2);
            sem_oric(ibin) = std(catchDirectionDeg(:,indc),[],2)./sqrt(length(indc));
            FAs_ori(ibin) = sum(FAIx(indc),2);
            CRs_ori(ibin) = sum(CRIx(indc),2);
%             rct_oric(ibin) = mean(catchReact(FAIx(indc)),2);
%             rctsem_oric(ibin) = std(catchReact(FAIx(indc)),[],2)./length(indc);
        end
    end
    [HR_ori, ci95_HR_ori] = binofit(Hits_ori, Misses_ori + Hits_ori);
    [FR_ori, ci95_FR_ori] = binofit(FAs_ori, CRs_ori + FAs_ori);
    n_ori = [Hits_ori+Misses_ori; FAs_ori+CRs_ori];
    
    
    
    
    figure; 
    errorbarxy(avg_ori, HR_ori, sem_ori, sem_ori, HR_ori - ci95_HR_ori(:,1)', ci95_HR_ori(:,2)' - HR_ori, {'ok', 'k', 'k'});
    hold on
    errorbarxy(avg_oric, FR_ori, sem_oric, sem_oric, FR_ori - ci95_FR_ori(:,1)', ci95_FR_ori(:,2)' - FR_ori, {'oc', 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['i' num2str(mouse_name) '- Visual and visual catch trials'])
    xlabel('Orientation change (deg)')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([10 100])
    legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);

    
    
    %% plot visual trials and visual catch hit rate by pulse number
    tCyclesOn = cell2mat(input.tCyclesOn);
catchCyclesOn = cell2mat(input.catchCyclesOn);
cycles = unique(tCyclesOn);

cycAvg_ori = zeros(length(cycles),max(ori_bins,[],2));
cycSem_ori = zeros(length(cycles),max(ori_bins,[],2));
cycHits_ori = zeros(length(cycles),max(ori_bins,[],2));
cycMisses_ori = zeros(length(cycles),max(ori_bins,[],2));

cycFAs_ori = zeros(length(cycles),max(ori_bins,[],2));
cycCRs_ori = zeros(length(cycles),max(ori_bins,[],2));
cycSem_oric = zeros(length(cycles),max(ori_bins,[],2));
cycAvg_oric = zeros(length(cycles),max(ori_bins,[],2));

cycHR_ori = zeros(length(cycles),max(ori_bins,[],2));
ci95_cycHR_ori = zeros(max(ori_bins,[],2),2,length(cycles));
cycFR_ori = zeros(length(cycles),max(ori_bins,[],2));
ci95_cycFR_ori = zeros(max(ori_bins,[],2),2,length(cycles));

for icyc = 1:length(cycles);
for ibin = 1:max(ori_bins,[],2)
        ind = intersect(find(tCyclesOn == cycles(icyc)),find(bin_ori == ibin));
        cycAvg_ori(icyc,ibin) = mean(targetDirectionDeg(:,ind),2);
        cycSem_ori(icyc,ibin) = std(targetDirectionDeg,[],2)./sqrt(length(ind));
        cycHits_ori(icyc,ibin) = sum(successIx(ind),2);
        cycMisses_ori(icyc,ibin) = sum(missedIx(ind),2);

        indc = intersect(find(catchCyclesOn == cycles(icyc)),find(bin_oric == ibin));
        cycAvg_oric(icyc,ibin) = mean(catchDirectionDeg(:,indc),2);
        cycSem_oric(icyc,ibin) = std(catchDirectionDeg(:,indc),[],2)./sqrt(length(indc));
        cycFAs_ori(icyc,ibin) = sum(FAIx(indc),2);
        cycCRs_ori(icyc,ibin) = sum(CRIx(indc),2);
end
[cycHR_ori(icyc,:), ci95_cycHR_ori(:,:,icyc)] = binofit(cycHits_ori(icyc,:),(cycMisses_ori(icyc,:) + cycHits_ori(icyc,:)));
[cycFR_ori(icyc,:), ci95_cycFR_ori(:,:,icyc)] = binofit(cycFAs_ori(icyc,:), cycCRs_ori(icyc,:) + cycFAs_ori(icyc,:));
end

ci95_cycHR_ori_U = squeeze(ci95_cycHR_ori(:,2,:))' - cycHR_ori;
ci95_cycHR_ori_L = cycHR_ori - squeeze(ci95_cycHR_ori(:,1,:))';
ci95_cycFR_ori_U = squeeze(ci95_cycFR_ori(:,2,:))' - cycFR_ori;
ci95_cycFR_ori_L = cycFR_ori - squeeze(ci95_cycFR_ori(:,1,:))';


cycN_HR_ori = [cycHits_ori+cycMisses_ori];
cycN_FA_ori = [cycFAs_ori+cycCRs_ori];

catchDegInd = find(~isnan(FR_ori));
nTrInd_FA = cycN_FA_ori(:,catchDegInd)>5;
nTrInd_HR = cycN_HR_ori(:,catchDegInd)>5;

% figure;
% subplot(1,2,1)
% errorbar(repmat(cycles(nTrInd_HR(:,1)),[length(catchDegInd(1)),1])', cycHR_ori((nTrInd_HR(:,1)),catchDegInd(1)),ci95_cycHR_ori_L(nTrInd_HR(:,1),catchDegInd(1)), ci95_cycHR_ori_U(nTrInd_HR(:,1),catchDegInd(1)),'ko-');
% hold on
% errorbar(repmat(cycles(nTrInd_FA(:,1)),[length(catchDegInd(1)),1])', cycFR_ori((nTrInd_FA(:,1)),catchDegInd(1)),ci95_cycFR_ori_L(nTrInd_FA(:,1),catchDegInd(1)), ci95_cycFR_ori_U(nTrInd_FA(:,1),catchDegInd(1)),'co-');
% ylim([0 1])
% xlim([0 11])
% title(['hard trials; ' num2str(avg_ori(:,catchDegInd(1))) '/' num2str(avg_oric(:,catchDegInd(1))) ' avg deg-valid/invalid'])
% xlabel('Cycle #')
% ylabel('Hit rate')
% legend({'valid'; 'invalid'},'Location','NorthWest')
% 
% subplot(1,2,2)
% errorbar(repmat(cycles(nTrInd_HR(:,2)),[length(catchDegInd(2)),1])', cycHR_ori((nTrInd_HR(:,2)),catchDegInd(2)),ci95_cycHR_ori_L(nTrInd_HR(:,2),catchDegInd(2)), ci95_cycHR_ori_U(nTrInd_HR(:,2),catchDegInd(2)),'ko-');
% hold on
% errorbar(repmat(cycles(nTrInd_FA(:,2)),[length(catchDegInd(2)),1])', cycFR_ori((nTrInd_FA(:,2)),catchDegInd(2)),ci95_cycFR_ori_L(nTrInd_FA(:,2),catchDegInd(2)), ci95_cycFR_ori_U(nTrInd_FA(:,2),catchDegInd(2)),'co-');
% ylim([0 1])
% xlim([0 11])
% title(['easy trials; ' num2str(avg_ori(:,catchDegInd(2))) '/' num2str(avg_oric(:,catchDegInd(2))) ' avg deg-valid/invalid'])
% xlabel('Cycle #')
% ylabel('Hit rate')
% 
% suptitle([num2str(mouse_name) ' - Hit rate as a function of trial length'])
% print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_HRandFRbyPulseN.pdf']),'-dpdf')
    %%
    figure(VFig);
    ind_catch = find(n_ori(2,:)>5);
    HitsV = zeros(1,length(ind_catch));
    MissesV = zeros(1,length(ind_catch));
    FAsV = zeros(1,length(ind_catch));
    CRsV = zeros(1,length(ind_catch));
    if ~isempty(ind_catch)
    for ibin = 1:length(ind_catch)
       tarC = unique(catchDirectionDeg(find(bin_oric==(ind_catch(ibin)))));
       indT = find(ismember(targetDirectionDeg,tarC));
       indC = find(ismember(catchDirectionDeg,tarC));
       HitsV(ibin) = sum(successIx(indT),2);
       MissesV(ibin) = sum(missedIx(indT),2);
       FAsV(ibin) = sum(FAIx(indC),2);
       CRsV(ibin) = sum(CRIx(indC),2); 
    end
    [HR_V, HR_ci95_V] = binofit(HitsV, HitsV+MissesV);
    [FR_V, FR_ci95_V] = binofit(FAsV, FAsV+CRsV);
    errorbarxy(HR_V, FR_V, HR_V-HR_ci95_V(:,1)', HR_ci95_V(:,2)'-HR_V, FR_V-FR_ci95_V(:,1)', FR_ci95_V(:,2)'-FR_V, {['o' av(ms_ind).col_str], av(ms_ind).col_str, av(ms_ind).col_str})
    hold on
    xlim([0 1])
    ylim([0 1])
    plot([0:0.1:1],[0:0.1:1],'--k')
    hold on
    xlabel('Hit rate- valid cue')
    ylabel('Hit rate- invalid cue')
    title('Visual trials')
    end

end


figure(mice_hr_fig)
subplot(1,2,1)
hold on
h = errorbar([1 2],mean(easy_all,2),std(easy_all,[],2)./sqrt(size(easy_all,2)),'ko-');
h.MarkerFaceColor = 'k';
h.LineWidth = 2;
ylim([0 1])
xlim([0.75 2.25])
mice_hr_fig.CurrentAxes.XTick = [1 2];
mice_hr_fig.CurrentAxes.XTickLabel = {'valid';'invalid'};
title('easy')
ylabel('HR')


subplot(1,2,2)
hold on
h = errorbar([1 2],mean(hard_all,2),std(hard_all,[],2)./sqrt(size(hard_all,2)),'ko-');
h.MarkerFaceColor = 'k';
h.LineWidth = 2;
ylim([0 1])
xlim([0.75 2.25])
mice_hr_fig.CurrentAxes.XTick = [1 2];
mice_hr_fig.CurrentAxes.XTickLabel = {'valid';'invalid'};
title('hard')
ylabel('HR')

print(fullfile(rc.fitOutputSummary, [date '_HRplotAllExpts_fsavImg.pdf']),'-dpdf')


figure(VFig);
print(fullfile(rc.fitOutputSummary, [date '_VisualCatchSummary_fsavImg.pdf']),'-dpdf')
figure(AFig);
print(fullfile(rc.fitOutputSummary, [date '_AuditoryCatchSummary_fsavImg.pdf']),'-dpdf')


%% plot HR summary across mice
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
set(0,'DefaultaxesFontSize', 16)
%minimum number of trials in order to plot
minTrNms = 1;
minTrNall = 1;
% create universal ori and amp bins
tarDeg_all = cat(2,cell2mat(tarDeg));
cDeg_all = cat(2,cell2mat(cDeg));
% tarAmp_all = cat(2,cell2mat(tarAmp));
% cAmp_all = cat(2,cell2mat(cAmp));
sIx = cat(2,cell2mat(sIx_all));
mIx = cat(2,cell2mat(mIx_all));
faIx = cat(2,cell2mat(faIx_all));
crIx = cat(2,cell2mat(crIx_all));
eIx = cat(2,cell2mat(eIx_all));
isCatchFACR = bsxfun(@plus,faIx,crIx);
cycTimeMs_all = cat(2,cell2mat(cycTimeMs_ms));
% % tarRct_all = cat(2,cell2mat(tarRct));
% % cRct_all = cat(2,cell2mat(cRct));
% % eRct_all = cat(2,cell2mat(eRct));
% prevTrVis_all = cat(2,cell2mat(prevTrVis));
% prevTrAud_all = cat(2,cell2mat(prevTrAud));

reducedOriC_edges = [1 24 100];
catch_ori_edges = unique(cDeg_all(~isnan(cDeg_all)));
% catch_amp_edges = unique(cAmp_all(~isnan(cAmp_all)));

[h_ori, bin_ori_all] = histc(tarDeg_all, ori_edges);
[h_oric, bin_oric_all] = histc(cDeg_all, ori_edges);
% [h_amp, bin_amp_all] = histc(tarAmp_all, amp_edges);
% [h_ampc, bin_ampc_all] = histc(cAmp_all, amp_edges);


nTrials = cell2mat(cellfun(@length,tarDeg,'Unif',false));
cumTrials = cumsum(nTrials);

ori_bins = unique(bin_ori_all);
% amp_bins = unique(bin_amp_all);
%% for each mouse with universal bins
% get hit rate curves across all mice with universal ori and amp bins
for imouse = 1:nMice
    
    if imouse == 1
        trialInd = 1:nTrials(imouse);
    else
        trialInd = cumTrials(imouse-1)+1:cumTrials(imouse);
    end
    %pre-allocate
    Hits_ori_all{imouse} = zeros(1,max(ori_bins,[],2));
    Misses_ori_all{imouse} = zeros(1,max(ori_bins,[],2));
    FAs_ori_all{imouse} = zeros(1,max(ori_bins,[],2));
    CRs_ori_all{imouse} = zeros(1,max(ori_bins,[],2));
    
%     Hits_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
%     Misses_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
%     FAs_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
%     CRs_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
%     
    n_ori_all{imouse} = zeros(2,max(ori_bins,[],2));
%     n_amp_all{imouse} = zeros(2,max(amp_bins,[],2));
    
    rct_ori_each{imouse} = NaN(1,max(ori_bins,[],2));
    rctsem_ori_each{imouse} = NaN(1,max(ori_bins,[],2));
    rct_oric_each{imouse} = NaN(1,max(ori_bins,[],2));
    rctsem_oric_each{imouse} = NaN(1,max(ori_bins,[],2));
%     rct_amp_each{imouse} = NaN(1,max(amp_bins,[],2));
%     rctsem_amp_each{imouse} = NaN(1,max(amp_bins,[],2));
%     rct_ampc_each{imouse} = NaN(1,max(amp_bins,[],2));
%     rctsem_ampc_each{imouse} = NaN(1,max(amp_bins,[],2));
    
    for ibin = 1:max(ori_bins,[],2)
        ind_v = find(bin_ori_all(trialInd) == ibin);
        if length(ind_v)>minTrNms
            Hits_ori_all{imouse}(ibin) = sum(sIx_all{imouse}(ind_v),2);
            Misses_ori_all{imouse}(ibin) = sum(mIx_all{imouse}(ind_v),2);
% %             rct_ori_each{imouse}(ibin) = mean(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind_v)),2);
% %             rctsem_ori_each{imouse}(ibin) = std(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind_v)),[],2)./sqrt(length(ind_v));
        end
        ind_av = find(bin_oric_all(trialInd) == ibin);
        if length(ind_av)>minTrNms
            FAs_ori_all{imouse}(ibin) = sum(faIx_all{imouse}(ind_av),2);
            CRs_ori_all{imouse}(ibin) = sum(crIx_all{imouse}(ind_av),2);
% %             rct_oric_each{imouse}(ibin) = mean(cRct{imouse}(intersect(find(faIx_all{imouse}),ind_av)),2);
% %             rctsem_oric_each{imouse}(ibin) = std(cRct{imouse}(intersect(find(faIx_all{imouse}),ind_av)),[],2)./sqrt(length(ind_av));
        end
    end
    [HR_ori_each{imouse}, ci95_HR_ori_each{imouse}] = binofit(Hits_ori_all{imouse}, Misses_ori_all{imouse} + Hits_ori_all{imouse});
    [FR_ori_each{imouse}, ci95_FR_ori_each{imouse}] = binofit(FAs_ori_all{imouse}, CRs_ori_all{imouse} + FAs_ori_all{imouse});
    n_ori_all{imouse} = [Hits_ori_all{imouse}+Misses_ori_all{imouse}; FAs_ori_all{imouse}+CRs_ori_all{imouse}];
    
% 
%     for ibin = 1:max(amp_bins,[],2)
%         ind_v = find(bin_amp_all(trialInd) == ibin);
%         if length(ind_v)>minTrNms
%             Hits_amp_all{imouse}(ibin) = sum(sIx_all{imouse}(ind_v),2);
%             Misses_amp_all{imouse}(ibin) = sum(mIx_all{imouse}(ind_v),2);
% % %             rct_amp_each{imouse}(ibin) = mean(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind_v)),2);
% % %             rctsem_amp_each{imouse}(ibin) = std(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind_v)),[],2)./sqrt(length(ind_v));
%         end
%         ind_av = find(bin_ampc_all(trialInd) == ibin);
%         if length(ind_av)>minTrNms
%             FAs_amp_all{imouse}(ibin) = sum(faIx_all{imouse}(ind_av),2);
%             CRs_amp_all{imouse}(ibin) = sum(crIx_all{imouse}(ind_av),2);
%             rct_ampc_each{imouse}(ibin) = mean(cRct{imouse}(intersect(find(faIx_all{imouse}),ind_av)),2);
%             rctsem_ampc_each{imouse}(ibin) = std(cRct{imouse}(intersect(find(faIx_all{imouse}),ind_av)),[],2)./sqrt(length(ind_av));
%         end
%     end
%     [HR_amp_each{imouse}, ci95_HR_amp_each{imouse}] = binofit(Hits_amp_all{imouse}, Misses_amp_all{imouse} + Hits_amp_all{imouse});
%     [FR_amp_each{imouse}, ci95_FR_amp_each{imouse}] = binofit(FAs_amp_all{imouse}, CRs_amp_all{imouse} + FAs_amp_all{imouse});
%     n_amp_all{imouse} = [Hits_amp_all{imouse}+Misses_amp_all{imouse}; FAs_amp_all{imouse}+CRs_amp_all{imouse}];
end

% visMiceInd = intersect(find(~isnan(nanmean(cell2mat(rct_ori_each'),2))),find(~isnan(nanmean(cell2mat(rct_oric_each'),2))));
% audMiceInd = intersect(find(~isnan(nanmean(cell2mat(rct_amp_each'),2))),find(~isnan(nanmean(cell2mat(rct_ampc_each'),2))));


%% use bins across all mice
% tarDeg_all = cat(2,cell2mat(tarDeg(visMiceInd)));
% cDeg_all = cat(2,cell2mat(cDeg(visMiceInd)));
% sIx = cat(2,cell2mat(sIx_all(visMiceInd)));
% mIx = cat(2,cell2mat(mIx_all(visMiceInd)));
% faIx = cat(2,cell2mat(faIx_all(visMiceInd)));
% crIx = cat(2,cell2mat(crIx_all(visMiceInd)));
% eIx = cat(2,cell2mat(eIx_all(visMiceInd)));
% % tarRct_all = cat(2,cell2mat(tarRct(visMiceInd)));
% % cRct_all = cat(2,cell2mat(cRct(visMiceInd)));
% % eRct_all = cat(2,cell2mat(eRct(visMiceInd)));
% prevTrVis_all = cat(2,cell2mat(prevTrVis(visMiceInd)));
% prevTrAud_all = cat(2,cell2mat(prevTrAud(visMiceInd)));

% tarAmp_all = cat(2,cell2mat(tarAmp(audMiceInd)));
% cAmp_all = cat(2,cell2mat(cAmp(audMiceInd)));
% sIx_A = cat(2,cell2mat(sIx_all(audMiceInd)));
% mIx_A = cat(2,cell2mat(mIx_all(audMiceInd)));
% faIx_A = cat(2,cell2mat(faIx_all(audMiceInd)));
% crIx_A = cat(2,cell2mat(crIx_all(audMiceInd)));
% eIx_A = cat(2,cell2mat(eIx_all(audMiceInd)));
% tarRct_A_all = cat(2,cell2mat(tarRct(audMiceInd)));
% cRct_A_all = cat(2,cell2mat(cRct(audMiceInd)));
% eRct_A_all = cat(2,cell2mat(eRct(audMiceInd)));
% % prevTrVis_A_all = cat(2,cell2mat(prevTrVis(audMiceInd)));
% % prevTrAud_A_all = cat(2,cell2mat(prevTrAud(audMiceInd)));

[h_ori, bin_ori_all] = histc(tarDeg_all, ori_edges);
[h_oric, bin_oric_all] = histc(cDeg_all, ori_edges);
% [h_amp, bin_amp_all] = histc(tarAmp_all, amp_edges);
% [h_ampc, bin_ampc_all] = histc(cAmp_all, amp_edges);

avg_ori_all = zeros(1,max(ori_bins,[],2));
sem_ori_all = zeros(1,max(ori_bins,[],2));    
% avg_amp_all = zeros(1,max(amp_bins,[],2));
% sem_amp_all = zeros(1,max(amp_bins,[],2));
avg_oric_all = zeros(1,max(ori_bins,[],2));
sem_oric_all = zeros(1,max(ori_bins,[],2));    
% avg_ampc_all = zeros(1,max(amp_bins,[],2));
% sem_ampc_all = zeros(1,max(amp_bins,[],2));
% % rct_ori_all = NaN(1,max(ori_bins,[],2));
% % rctsem_ori_all = NaN(1,max(ori_bins,[],2));
% % rct_oric_all = NaN(1,max(ori_bins,[],2));
% % rctsem_oric_all = NaN(1,max(ori_bins,[],2));
% % % rct_amp_all = NaN(1,max(amp_bins,[],2));
% % % rctsem_amp_all = NaN(1,max(amp_bins,[],2));
% % % rct_ampc_all = NaN(1,max(amp_bins,[],2));
% % % rctsem_ampc_all = NaN(1,max(amp_bins,[],2));
ori_bin_indices = cell(1,max(ori_bins,[],2));
oric_bin_indices = cell(1,max(ori_bins,[],2));
% amp_bin_indices = cell(1,max(ori_bins,[],2));
% ampc_bin_indices = cell(1,max(ori_bins,[],2));

% % variables for dependence on previous trial
% avg_ori_prevVis = zeros(1,max(ori_bins,[],2));
% sem_ori_prevVis = zeros(1,max(ori_bins,[],2)); 
% avg_oric_prevVis = zeros(1,max(ori_bins,[],2));
% sem_oric_prevVis = zeros(1,max(ori_bins,[],2));  
% avg_ori_prevAud = zeros(1,max(ori_bins,[],2));
% sem_ori_prevAud = zeros(1,max(ori_bins,[],2)); 
% avg_oric_prevAud = zeros(1,max(ori_bins,[],2));
% sem_oric_prevAud = zeros(1,max(ori_bins,[],2));  

for ibin = 1:max(ori_bins,[],2)
    ind_v = find(bin_ori_all == ibin);
    ori_bin_indices{ibin} = ind_v;
    avg_ori_all(ibin) = mean(tarDeg_all(ind_v),2);
    sem_ori_all(ibin) = std(tarDeg_all(ind_v),[],2)./sqrt(length(ind_v));
% %     rct_ori_all(ibin) = mean(tarRct_all(intersect(find(sIx),ind_v)),2);
% %     rctsem_ori_all(ibin) = std(tarRct_all(intersect(find(sIx),ind_v)),[],2)./sqrt(length(ind_v));
    ind_av = find(bin_oric_all == ibin);
    oric_bin_indices{ibin} = ind_av;
    avg_oric_all(ibin) = mean(cDeg_all(ind_av),2);
    sem_oric_all(ibin) = std(cDeg_all(ind_av),[],2)./sqrt(length(ind_av));
% %     rct_oric_all(ibin) = mean(cRct_all(intersect(find(faIx),ind_av)),2);
% %     rctsem_oric_all(ibin) = std(cRct_all(intersect(find(faIx),ind_av)),[],2)./sqrt(length(ind_av));
    
    % dependence on previous trial
%     pv_ind = intersect(ind_v,find(prevTrVis_all));
%     pa_ind = intersect(ind_v,find(prevTrAud_all));
%     avg_ori_prevVis(ibin) = mean(tarDeg_all(pv_ind),2);
%     sem_ori_prevVis(ibin) = std(tarDeg_all(pv_ind),[],2)./sqrt(length(pv_ind));
%     avg_ori_prevAud(ibin) = mean(tarDeg_all(pa_ind),2);
%     sem_ori_prevAud(ibin) = std(tarDeg_all(pa_ind),[],2)./sqrt(length(pa_ind));
%     pv_indc = intersect(ind_av,find(prevTrVis_all));
%     pa_indc = intersect(ind_av,find(prevTrAud_all));
%     avg_oric_prevVis(ibin) = mean(cDeg_all(pv_indc),2);
%     sem_oric_prevVis(ibin) = std(cDeg_all(pa_indc),[],2)./sqrt(length(pa_indc));
%     avg_oric_prevAud(ibin) = mean(cDeg_all(pa_indc),2);
%     sem_oric_prevAud(ibin) = std(cDeg_all(pa_indc),[],2)./sqrt(length(pa_indc));
end
% 
% for ibin = 1:max(amp_bins,[],2);
%     ind_v = find(bin_amp_all == ibin);
%     amp_bin_indices{ibin} = ind_v;
%     avg_amp_all(ibin) = mean(tarAmp_all(ind_v),2);
%     sem_amp_all(ibin) = std(tarAmp_all(ind_v),[],2)./sqrt(length(ind_v));
%     rct_amp_all(ibin) = mean(tarRct_A_all(intersect(find(sIx_A),ind_v)),2);
%     rctsem_amp_all(ibin) = std(tarRct_A_all(intersect(find(sIx_A),ind_v)),[],2)./sqrt(length(ind_v));
%     ind_av = find(bin_ampc_all == ibin);
%     ampc_bin_indices{ibin} = ind_av;
%     avg_ampc_all(ibin) = mean(cAmp_all(ind_av),2);
%     sem_ampc_all(ibin) = std(cAmp_all(ind_av),[],2)./sqrt(length(ind_av));
%     rct_ampc_all(ibin) = mean(cRct_A_all(intersect(find(faIx_A),ind_av)),2);
%     rctsem_ampc_all(ibin) =  std(cRct_A_all(intersect(find(faIx_A),ind_av)),[],2)./sqrt(length(ind_av));
% end


hitsAll_V = zeros(1,max(ori_bins,[],2));
missAll_V = zeros(1,max(ori_bins,[],2));
faAll_V = zeros(1,max(ori_bins,[],2));
crAll_V = zeros(1,max(ori_bins,[],2));
n_ori = zeros(1,max(ori_bins,[],2));
% 
% hitsAll_V_prevV = zeros(1,max(ori_bins,[],2));
% missAll_V_prevV = zeros(1,max(ori_bins,[],2));
% faAll_V_prevV = zeros(1,max(ori_bins,[],2));
% crAll_V_prevV = zeros(1,max(ori_bins,[],2));
% n_ori_prevV = zeros(1,max(ori_bins,[],2));

% hitsAll_V_prevA = zeros(1,max(ori_bins,[],2));
% missAll_V_prevA = zeros(1,max(ori_bins,[],2));
% faAll_V_prevA = zeros(1,max(ori_bins,[],2));
% crAll_V_prevA = zeros(1,max(ori_bins,[],2));
% n_ori_prevA = zeros(1,max(ori_bins,[],2));

for ibin = 1:max(ori_bins,[],2)
    ind_v = find(bin_ori_all == ibin);
%     pv_ind = intersect(ind_v,find(prevTrVis_all));
%     pa_ind = intersect(ind_v,find(prevTrAud_all));
    if sum(sIx(ind_v)) + sum(mIx(ind_v)) > minTrNall
    hitsAll_V(ibin) = sum(sIx(ind_v));
    missAll_V(ibin)= sum(mIx(ind_v));
    
%     hitsAll_V_prevV(ibin) = sum(sIx(pv_ind));
%     missAll_V_prevV(ibin) = sum(mIx(pv_ind));
%     
%     hitsAll_V_prevA(ibin) = sum(sIx(pa_ind));
%     missAll_V_prevA(ibin) = sum(mIx(pa_ind));
%     
    end
    ind_av = find(bin_oric_all == ibin);
%     pv_indc = intersect(ind_av,find(prevTrVis_all));
%     pa_indc = intersect(ind_av,find(prevTrAud_all));
    if sum(faIx(ind_av)) + sum(crIx(ind_av)) > minTrNall
    faAll_V(ibin) = sum(faIx(ind_av));
    crAll_V(ibin) = sum(crIx(ind_av)); 
% 
%     faAll_V_prevV(ibin) = sum(faIx(pv_indc));
%     crAll_V_prevV(ibin) = sum(crIx(pv_indc)); 
%     
%     faAll_V_prevA(ibin) = sum(faIx(pa_indc));
%     crAll_V_prevA(ibin) = sum(crIx(pa_indc)); 
%     
    end
    n_ori(ibin) = hitsAll_V(ibin)+missAll_V(ibin);
    n_oric(ibin) = faAll_V(ibin)+crAll_V(ibin);
%     n_ori_prevV(ibin) = hitsAll_V_prevV(ibin)+missAll_V_prevV(ibin);
%     n_oric_prevV(ibin) = faAll_V_prevV(ibin)+crAll_V_prevV(ibin);
%     n_ori_prevA(ibin) = hitsAll_V_prevA(ibin)+missAll_V_prevA(ibin);
%     n_oric_prevA(ibin) = faAll_V_prevA(ibin)+crAll_V_prevA(ibin);
end

[HR_ori_all, ci95_HR_ori_all] = binofit(hitsAll_V, missAll_V + hitsAll_V);
[FR_ori_all, ci95_FR_ori_all] = binofit(faAll_V, faAll_V + crAll_V);

% [HR_ori_prevV, ci95_HR_ori_prevV] = binofit(hitsAll_V_prevV, missAll_V_prevV + hitsAll_V_prevV);
% [FR_ori_prevV, ci95_FR_ori_prevV] = binofit(faAll_V_prevV, faAll_V_prevV + crAll_V_prevV);
% 
% [HR_ori_prevA, ci95_HR_ori_prevA] = binofit(hitsAll_V_prevA, missAll_V_prevA + hitsAll_V_prevA);
% [FR_ori_prevA, ci95_FR_ori_prevA] = binofit(faAll_V_prevA, faAll_V_prevA + crAll_V_prevA);

% n_amp = zeros(1,max(amp_bins,[],2));
% n_ampc = zeros(1,max(amp_bins,[],2));
% hitsAll_A = zeros(1,max(amp_bins,[],2));
% missAll_A= zeros(1,max(amp_bins,[],2));
% faAll_A = zeros(1,max(amp_bins,[],2));
% crAll_A= zeros(1,max(amp_bins,[],2));
% for ibin = 1:max(amp_bins,[],2)
%     ind_v = find(bin_amp_all == ibin);
%     if  sum(sIx_A(ind_v)) + sum(mIx_A(ind_v)) >minTrNall
%     hitsAll_A(ibin) = sum(sIx_A(ind_v));
%     missAll_A(ibin)= sum(mIx_A(ind_v));
%     end
%     ind_av = find(bin_ampc_all == ibin);
%     if  sum(faIx_A(ind_av)) + sum(crIx_A(ind_av)) >minTrNall
%     faAll_A(ibin) = sum(faIx_A(ind_av));
%     crAll_A(ibin) = sum(crIx_A(ind_av));
%     end
%     n_amp(ibin) = hitsAll_A(ibin)+missAll_A(ibin);
%     n_ampc(ibin) = faAll_A(ibin)+crAll_A(ibin);
% end
% 
% [HR_amp_all, ci95_HR_amp_all] = binofit(hitsAll_A, missAll_A + hitsAll_A);
% [FR_amp_all, ci95_FR_amp_all] = binofit(faAll_A, faAll_A + crAll_A);


%% plot each mouse curve and avg across mice with universal bins
uniHRFig = figure;
uniRTFig = figure;


% xAud = [min(avg_amp_all)/2 1];
xVis = [min(avg_ori_all)/2 100];
yHR = [0 1];
yRT = [0 550];

xV1 = [floor(log10(xVis(1))) ceil(log10(xVis(end)))];
xVTick = 10.^(xV1(1):1:xV1(2));
XVTick = xVTick(xVTick>=xVis(1) & xVTick <=xVis(end));
if length(xVTick) == 0
    xVTick = chop(xVis,2);
end
xVLabel = cellstr(num2str(xVTick(:)));
% 
% xA1 = [floor(log10(xAud(1))) ceil(log10(xAud(end)))];
% xATick = 10.^(xA1(1):1:xA1(2));
% XVTick = xATick(xATick>=xAud(1) & xATick <=xAud(end));
% if length(xATick) == 0
%     xATick = chop(xAud,2);
% end
% xALabel = cellstr(num2str(xATick(:)));

%visual trials
for imouse = 1:nMice
    figure(uniHRFig)
    subplot(2,1,1)
    hold on
    colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
    hrInd = ~isnan(HR_ori_each{imouse});
    frInd = ~isnan(FR_ori_each{imouse});
    plot(avg_ori_all(hrInd),HR_ori_each{imouse}(hrInd),'ko-')
    hold on
    plot(avg_ori_all(frInd),FR_ori_each{imouse}(frInd),'co-')
    hold on
%     errorbarxy(avg_ori_all, HR_ori_each{imouse}, sem_ori_all, sem_ori_all, HR_ori_each{imouse} - ci95_HR_ori_each{imouse}(:,1)', ci95_HR_ori_each{imouse}(:,2)' - HR_ori_each{imouse}, {['o' 'k'], 'k', 'k'});
%     hold on
%     errorbarxy(avg_oric_all, FR_ori_each{imouse}, sem_oric_all, sem_oric_all, FR_ori_each{imouse} - ci95_FR_ori_each{imouse}(:,1)', ci95_FR_ori_each{imouse}(:,2)' - FR_ori_each{imouse}, {['o' 'c'], 'c', 'c'});
%     set(gca, 'xscale', 'log');
    title(['visual trials'])
    xlabel('Orientation change (deg)')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([0 100])
%     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
    
end
figure(uniHRFig)
subplot(2,1,1)
set(gca, 'YLim', yHR, ...
       'XLim', xVis, ...
       'XScale', 'log', ...
       'XTick', xVTick, ...
       'XTickLabel', xVLabel);
% xlim(xVis)
% ylim(yHR)
% set(gca,'XScale','log')
hold on
errorbarxy(avg_ori_all, HR_ori_all, sem_ori_all, sem_ori_all, HR_ori_all - ci95_HR_ori_all(:,1)', ci95_HR_ori_all(:,2)' - HR_ori_all, {['o' 'k'], 'k', 'k'});
hold on
errorbarxy(avg_oric_all, FR_ori_all, sem_oric_all, sem_oric_all, FR_ori_all - ci95_FR_ori_all(:,1)', ci95_FR_ori_all(:,2)' - FR_ori_all, {['o' 'b'], 'b', 'b'});
hold on
hrInd = intersect(find(~isnan(HR_ori_all)),find(n_ori>minTrNall));
frInd = intersect(find(~isnan(FR_ori_all)),find(n_oric>minTrNall));
plot(avg_ori_all(hrInd),HR_ori_all(hrInd),'k','linewidth',3)
hold on
plot(avg_oric_all(frInd),FR_ori_all(frInd),'b','linewidth',3)
hold on
scatter(avg_ori_all(hrInd),HR_ori_all(hrInd),'ko','filled');
hold on
scatter(avg_oric_all(frInd),FR_ori_all(frInd),'bo','filled')
axis square

% %visual trials - RT
% for imouse = 1:nMice
%     figure(uniRTFig)
%     subplot(2,1,1)
%     hold on
%     colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
%     hrInd = ~isnan(rct_ori_each{imouse});
%     frInd = ~isnan(rct_oric_each{imouse});
%     plot(avg_ori_all(hrInd),rct_ori_each{imouse}(hrInd),'ko-')
%     hold on
%     plot(avg_oric_all(frInd),rct_oric_each{imouse}(frInd),'co-')
%     hold on
% %     errorbarxy(avg_ori_all, rct_ori_each{imouse}, sem_ori_all, sem_ori_all, rct_ori_each{imouse} - rctsem_ori_each{imouse}, rctsem_ori_each{imouse} - rct_ori_each{imouse}, {['o' 'k'], 'k', 'k'});
% %     hold on
% %     errorbarxy(avg_oric_all, rct_oric_each{imouse}, sem_oric_all, sem_oric_all, rct_oric_each{imouse} - rctsem_ori_each{imouse}, rctsem_oric_each{imouse} - rct_oric_each{imouse}, {['o' 'c'], 'c', 'c'});
% %     set(gca, 'xscale', 'log');
%     title(['visual trials'])
%     xlabel('Orientation change (deg)')
%     ylabel('RT')
% %     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
%     
% end
% figure(uniRTFig)
% subplot(2,1,1)
% set(gca, 'YLim', yRT, ...
%        'XLim', xVis, ...
%        'XScale', 'log', ...
%        'XTick', xVTick, ...
%        'XTickLabel', xVLabel);
% % xlim(xVis)
% % ylim(yRT)
% hold on
% errorbarxy(avg_ori_all, rct_ori_all, sem_ori_all, sem_ori_all, rct_ori_all - rctsem_ori_all, rctsem_ori_all - rct_ori_all, {['o' 'g'], 'g', 'g'});
% hold on
% errorbarxy(avg_oric_all, rct_oric_all, sem_oric_all, sem_oric_all, rct_oric_all - rctsem_oric_all, rctsem_oric_all - rct_oric_all, {['o' 'b'], 'b', 'b'});
% hold on
% hrInd = intersect(find(~isnan(rct_ori_all)),find(n_ori>minTrNall));
% frInd = intersect(find(~isnan(rct_oric_all)),find(n_oric>minTrNall));
% plot(avg_ori_all(hrInd),rct_ori_all(hrInd),'k','linewidth',3)
% hold on
% plot(avg_oric_all(frInd),rct_oric_all(frInd),'b','linewidth',3)
% hold on
% scatter(avg_ori_all(hrInd),rct_ori_all(hrInd),'go','filled');
% hold on
% scatter(avg_oric_all(frInd),rct_oric_all(frInd),'bo','filled');
% hold on
% axis square

% %auditory trials
% for imouse = 1:nMice
%     figure(uniHRFig)
%     subplot(2,1,2)
%     hold on
%     colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
%     hrInd = ~isnan(HR_amp_each{imouse});
%     frInd = ~isnan(FR_amp_each{imouse});
%     if length(hrInd) >= 2
%     plot(avg_amp_all(hrInd),HR_amp_each{imouse}(hrInd),'ko-')
%     hold on
%     end
%     if length(frInd) >=2
%     plot(avg_amp_all(frInd),FR_amp_each{imouse}(frInd),'co-')
%     hold on
%     end
% %     errorbarxy(avg_amp_all, HR_amp_each{imouse}, sem_amp_all, sem_amp_all, HR_amp_each{imouse} - ci95_HR_amp_each{imouse}(:,1)', ci95_HR_amp_each{imouse}(:,2)' - HR_amp_each{imouse}, {['o' 'k'], 'k', 'k'});
% %     hold on
% %     errorbarxy(avg_ampc_all, FR_amp_each{imouse}, sem_ampc_all, sem_ampc_all, FR_amp_each{imouse} - ci95_FR_amp_each{imouse}(:,1)', ci95_FR_amp_each{imouse}(:,2)' - FR_amp_each{imouse}, {['o' 'c'], 'c', 'c'});
% %     set(gca, 'xscale', 'log');
%     title(['auditory trials'])
%     xlabel('Target Sound Amplitude')
%     ylabel('Hit rate')
%     ylim([0 1])
%     xlim([0 1])
% %     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
%     
% end
% figure(uniHRFig)
% subplot(2,1,2)
% set(gca, 'YLim', yHR, ...
%        'XLim', xAud, ...
%        'XScale', 'log', ...
%        'XTick', xATick, ...
%        'XTickLabel', xALabel);
% % set(gca,'XTick',amp_edges)
% % xlim(xAud)
% % ylim(yHR)
% hold on
% errorbarxy(avg_amp_all, HR_amp_all, sem_amp_all, sem_amp_all, HR_amp_all - ci95_HR_amp_all(:,1)', ci95_HR_amp_all(:,2)' - HR_amp_all, {['o' 'k'], 'k', 'k'});
% hold on
% errorbarxy(avg_ampc_all, FR_amp_all, sem_ampc_all, sem_ampc_all, FR_amp_all - ci95_FR_amp_all(:,1)', ci95_FR_amp_all(:,2)' - FR_amp_all, {['o' 'b'], 'b', 'b'});
% hold on
% hrInd = intersect(find(~isnan(HR_amp_all)),find(n_amp>10));
% frInd = intersect(find(~isnan(HR_amp_all)),find(n_ampc>10));
% plot(avg_amp_all(hrInd),HR_amp_all(hrInd),'k','linewidth',3)
% hold on
% plot(avg_ampc_all(frInd),FR_amp_all(frInd),'b','linewidth',3)
% hold on
% scatter(avg_amp_all(hrInd),HR_amp_all(hrInd),'ko','filled');
% hold on
% scatter(avg_ampc_all(frInd),FR_amp_all(frInd),'bo','filled')
% axis square
% 
% %auditory trials - RT
% for imouse = 1:nMice
%     figure(uniRTFig)
%     subplot(2,1,2)
%     hold on
%     colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
%     hrInd = ~isnan(rct_amp_each{imouse});
%     frInd = ~isnan(rct_ampc_each{imouse});
%     if length(hrInd) >= 2
%     plot(avg_amp_all(hrInd),rct_amp_each{imouse}(hrInd),'ko--')
%     hold on
%     end
%     if length(frInd >=2)
%     plot(avg_amp_all(frInd),rct_amp_each{imouse}(frInd),'co--')
%     hold on
%     end
% %     errorbarxy(avg_ori_all, rct_ori_each{imouse}, sem_ori_all, sem_ori_all, rct_ori_each{imouse} - rctsem_ori_each{imouse}, rctsem_ori_each{imouse} - rct_ori_each{imouse}, {['o' 'k'], 'k', 'k'});
% %     hold on
% %     errorbarxy(avg_oric_all, rct_oric_each{imouse}, sem_oric_all, sem_oric_all, rct_oric_each{imouse} - rctsem_ori_each{imouse}, rctsem_oric_each{imouse} - rct_oric_each{imouse}, {['o' 'c'], 'c', 'c'});
% %     set(gca, 'xscale', 'log');
%     title(['auditory trials'])
%     xlabel('target sound amplitude')
%     ylabel('RT')
% %     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
%     
% end
% figure(uniRTFig)
% subplot(2,1,2)
% set(gca, 'YLim', yRT, ...
%        'XLim', xAud, ...
%        'XScale', 'log', ...
%        'XTick', xATick, ...
%        'XTickLabel', xALabel);
% % set(gca,'XTick',amp_edges)
% % xlim(xAud)
% % ylim(yRT)
% hold on
% errorbarxy(avg_amp_all, rct_amp_all, sem_amp_all, sem_amp_all, rct_amp_all - rctsem_amp_all, rctsem_amp_all - rct_amp_all, {['o' 'k'], 'k', 'k'});
% hold on
% errorbarxy(avg_ampc_all, rct_ampc_all, sem_ampc_all, sem_ampc_all, rct_ampc_all - rctsem_ampc_all, rctsem_ampc_all - rct_ampc_all, {['o' 'b'], 'b', 'b'});
% hold on
% hrInd = intersect(find(~isnan(rct_amp_all)),find(n_amp>minTrNall));
% frInd = intersect(find(~isnan(rct_ampc_all)),find(n_ampc>minTrNall));
% plot(avg_amp_all(hrInd),rct_amp_all(hrInd),'k--','linewidth',3)
% hold on
% plot(avg_ampc_all(frInd),rct_ampc_all(frInd),'b--','linewidth',3)
% hold on
% scatter(avg_amp_all(hrInd),rct_amp_all(hrInd),'ko','filled');
% hold on
% scatter(avg_ampc_all(frInd),rct_ampc_all(frInd),'bo','filled');
% hold on
% axis square

% %% dependence on previous trial figure
% set(0,'defaultfigurepaperorientation','landscape');
% set(0,'defaultfigurepapersize',[11 8.5]);
% set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
% set(0,'DefaultaxesFontSize', 16)
% uniHRFig_depPrevTr = figure;
% suptitle('visual trials - dep on prev trial')
% 
% hrInd = intersect(find(~isnan(HR_ori_all)),find(n_ori>minTrNall));
% frInd = intersect(find(~isnan(FR_ori_all)),find(n_oric>minTrNall));
% 
% %visual trials
% figure(uniHRFig_depPrevTr)
% subplot(2,3,1)
% hold on
% l1 = errorbarxy(avg_ori_all(hrInd), HR_ori_all(hrInd), sem_ori_all(hrInd), sem_ori_all(hrInd), HR_ori_all(hrInd) - ci95_HR_ori_all((hrInd),1)', ci95_HR_ori_all((hrInd),2)' - HR_ori_all(hrInd), {['o' 'k'], 'k', 'k'});
% l1.hMain.LineStyle = '-';
% l1.hMain.LineWidth = 3;
% l1.hMain.MarkerFaceColor = 'k';
% hold on
% l2 = errorbarxy(avg_oric_all(frInd), FR_ori_all(frInd), sem_oric_all(frInd), sem_oric_all(frInd), FR_ori_all(frInd) - ci95_FR_ori_all((frInd),1)', ci95_FR_ori_all((frInd),2)' - FR_ori_all(frInd), {['o' 'c'], 'c', 'c'});
% l2.hMain.LineStyle = '-';
% l2.hMain.LineWidth = 3;
% l2.hMain.MarkerFaceColor = 'c';
% hold on
% axis square
% title('all trials')
% legend([l1.hMain l2.hMain],'val','inv','Location','northeastoutside')
% 
% figure(uniHRFig_depPrevTr)
% subplot(2,3,2)
% hold on
% l1 = errorbarxy(avg_ori_prevVis(hrInd), HR_ori_prevV(hrInd), sem_ori_prevVis(hrInd), sem_ori_prevVis(hrInd), HR_ori_prevV(hrInd) - ci95_HR_ori_prevV((hrInd),1)', ci95_HR_ori_prevV((hrInd),2)' - HR_ori_prevV(hrInd), {['o' 'k'], 'k', 'k'});
% l1.hMain.LineStyle = '-';
% l1.hMain.LineWidth = 3;
% l1.hMain.MarkerFaceColor = 'k';
% hold on
% l2 = errorbarxy(avg_oric_prevVis(frInd), FR_ori_prevV(frInd), sem_oric_prevVis(frInd), sem_oric_prevVis(frInd), FR_ori_prevV(frInd) - ci95_FR_ori_prevV((frInd),1)', ci95_FR_ori_prevV((frInd),2)' - FR_ori_prevV(frInd), {['o' 'c'], 'c', 'c'});
% l2.hMain.LineStyle = '-';
% l2.hMain.LineWidth = 3;
% l2.hMain.MarkerFaceColor = 'c';
% hold on
% axis square
% title('previous trial = visual')
% legend([l1.hMain l2.hMain],'val','inv','Location','northeastoutside')
% 
% figure(uniHRFig_depPrevTr)
% subplot(2,3,3)
% hold on
% l1 = errorbarxy(avg_ori_prevAud(hrInd), HR_ori_prevA(hrInd), sem_ori_prevAud(hrInd), sem_ori_prevAud(hrInd), HR_ori_prevA(hrInd) - ci95_HR_ori_prevA((hrInd),1)', ci95_HR_ori_prevA((hrInd),2)' - HR_ori_prevA(hrInd), {['o' 'k'], 'k', 'k'});
% l1.hMain.LineStyle = '-';
% l1.hMain.LineWidth = 3;
% l1.hMain.MarkerFaceColor = 'k';
% hold on
% l2 = errorbarxy(avg_oric_prevAud(frInd), FR_ori_prevA(frInd), sem_oric_prevAud(frInd), sem_oric_prevAud(frInd), FR_ori_prevA(frInd) - ci95_FR_ori_prevA((frInd),1)', ci95_FR_ori_prevA((frInd),2)' - FR_ori_prevA(frInd), {['o' 'c'], 'c', 'c'});
% l2.hMain.LineStyle = '-';
% l2.hMain.LineWidth = 3;
% l2.hMain.MarkerFaceColor = 'c';
% hold on
% axis square
% title('previous trial = auditory')
% legend([l1.hMain l2.hMain],'val','inv','Location','northeastoutside')
% 
% figure(uniHRFig_depPrevTr)
% subplot(2,3,4)
% hold on
% l1 = errorbarxy(avg_ori_prevVis(hrInd), HR_ori_prevV(hrInd), sem_ori_prevVis(hrInd), sem_ori_prevVis(hrInd), HR_ori_prevV(hrInd) - ci95_HR_ori_prevV((hrInd),1)', ci95_HR_ori_prevV((hrInd),2)' - HR_ori_prevV(hrInd), {['o' 'k'], 'k', 'k'});
% l1.hMain.LineStyle = '-';
% l1.hMain.LineWidth = 3;
% l1.hMain.MarkerFaceColor = 'k';
% hold on
% l2 = errorbarxy(avg_ori_prevAud(hrInd), HR_ori_prevA(hrInd), sem_ori_prevAud(hrInd), sem_ori_prevAud(hrInd), HR_ori_prevA(hrInd) - ci95_HR_ori_prevA((hrInd),1)', ci95_HR_ori_prevA((hrInd),2)' - HR_ori_prevA(hrInd), {['o' 'k'], 'k', 'k'});
% l2.hMain.Color = [0.5 0.5 0.5];
% l2.hMain.LineStyle = '-';
% l2.hMain.LineWidth = 3;
% l2.hMain.MarkerFaceColor = [0.5 0.5 0.5];
% hold on
% axis square
% title('valid trials')
% legend([l1.hMain l2.hMain],'p vis','p aud','Location','northeastoutside')
% 
% figure(uniHRFig_depPrevTr)
% subplot(2,3,5)
% hold on
% l1 = errorbarxy(avg_oric_prevVis(frInd), FR_ori_prevV(frInd), sem_oric_prevVis(frInd), sem_oric_prevVis(frInd), FR_ori_prevV(frInd) - ci95_FR_ori_prevV((frInd),1)', ci95_FR_ori_prevV((frInd),2)' - FR_ori_prevV(frInd), {['o' 'c'], 'c', 'c'});
% l1.hMain.LineStyle = '-';
% l1.hMain.LineWidth = 3;
% l1.hMain.MarkerFaceColor = 'k';
% hold on
% l2 = errorbarxy(avg_oric_prevAud(frInd), FR_ori_prevA(frInd), sem_oric_prevAud(frInd), sem_oric_prevAud(frInd), FR_ori_prevA(frInd) - ci95_FR_ori_prevA((frInd),1)', ci95_FR_ori_prevA((frInd),2)' - FR_ori_prevA(frInd), {['o' 'c'], 'c', 'c'});
% l2.hMain.Color = [0 0.5 0.5];
% l2.hMain.LineStyle = '-';
% l2.hMain.LineWidth = 3;
% l2.hMain.MarkerFaceColor = [0 0.5 0.5];
% hold on
% axis square
% title('invalid trials')
% legend([l1.hMain l2.hMain],'p vis','p aud','Location','northeastoutside')
% 
% for iplot = 1:6
%     subplot(2,3,iplot)
%     
%     set(gca, 'YLim', yHR, ...
%            'XLim', xVis, ...
%            'XScale', 'log', ...
%            'XTick', xVTick, ...
%            'XTickLabel', xVLabel);
%        
%    xlabel('orientation')
%    ylabel('HR')
%    axis square
% end
% 
% print(fullfile(rc.fitOutputSummary, [date '_uniHR_depPrevTriald.pdf']),'-dpdf')
%% working memory test

% converge all react times for hits, invalid hits, earlies
visInd_all = find(tarDeg_all > 1);
cVisInd_all = find(cDeg_all > 1);
% audInd_all = find(tarAmp_all > 0);
% cAudInd_all = find(cAmp_all > 0);
% cRct_a_all = cat(2,cell2mat(cRct(audMiceInd)));
% eRct_a_all = cat(2,cell2mat(eRct(audMiceInd)));
% rct_a_h = mean(tarRct_all(intersect(audInd_all,find(cat(2,cell2mat(sIx_all(audMiceInd)))))));
% rct_a_h_sem = std(tarRct_all(intersect(audInd_all,find(cat(2,cell2mat(sIx_all(audMiceInd)))))))/sqrt(length(tarRct_all(intersect(audInd_all,find(cat(2,cell2mat(sIx_all(audMiceInd))))))));
% rct_a_fa = mean(cRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(faIx_all(audMiceInd)))))));
% rct_a_fa_sem = std(cRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(faIx_all(audMiceInd)))))))/sqrt(length(cRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(faIx_all(audMiceInd))))))));
% rct_a_e = nanmean(eRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(eIx_all(audMiceInd)))))));
% rct_a_e_sem = nanstd(eRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(eIx_all(audMiceInd)))))))/sqrt(sum(~isnan(eRct_a_all(intersect(audInd_all,find(cat(2,cell2mat(eIx_all(audMiceInd)))))))));
% rct_v_h = mean(tarRct_all(intersect(visInd_all,find(sIx))));
% rct_v_h_sem = std(tarRct_all(intersect(visInd_all,find(sIx))))/sqrt(length(tarRct_all(intersect(visInd_all,find(sIx)))));
% rct_v_fa = mean(cRct_all(intersect(cVisInd_all,find(faIx))));
% rct_v_fa_sem = std(cRct_all(intersect(cVisInd_all,find(faIx))))/sqrt(length(cRct_all(intersect(cVisInd_all,find(faIx)))));
% rct_v_e_sem = nanstd(eRct_all(intersect(visInd_all,find(eIx))))/sqrt(sum(~isnan(eRct_all(intersect(visInd_all,find(eIx))))));
% rct_v_e = nanmean(eRct_all(intersect(visInd_all,find(eIx))));


% rtAll = figure;
% subplot(2,1,1)
% errorbar(1,rct_v_h,rct_v_h_sem,'go');
% hold on
% errorbar(1,rct_v_fa,rct_v_fa_sem,'co');
% hold on
% errorbar(1,rct_v_e,rct_v_e_sem,'mo');
% xlabel('trial length (ms)')
% ylabel('RT')
% xlim([0 10])
% ylim([0 550])
% axis square
% title('visual trials - RT all')
% legend({'valid';'invalid';'earlies'})
% subplot(2,1,2)
% errorbar(1,rct_a_h,rct_a_h_sem,'ko');
% hold on
% errorbar(1,rct_a_fa,rct_a_fa_sem,'co');
% hold on
% errorbar(1,rct_a_e,rct_a_e_sem,'mo');
% xlabel('trial length (ms)')
% ylabel('RT')
% xlim([0 10])
% ylim([0 550])
% axis square
% title('auditory trials - RT all')
% legend({'valid';'invalid';'earlies'})

nTimes = 6;
nCyc_all =  double(cat(2,cell2mat(nCyc)));
trialTimeMs_all = cat(2,cell2mat(trialTimeMs));
trialTime_edges = [0:max(trialTimeMs_all)/nTimes:max(trialTimeMs_all)]+100;

[h_time, bin_time] = histc(trialTimeMs_all,trialTime_edges);


time_bins = unique(bin_time);
avg_trTime = zeros(1,length(time_bins));
hits_trTime = zeros(1,length(time_bins));
misses_trTime = zeros(1,length(time_bins));
earlies_trTime = zeros(1,length(time_bins));
CRbase_trTime = zeros(1,length(time_bins));
fa_trTime = zeros(1,length(time_bins));
cr_trTime = zeros(1,length(time_bins));
% rct_trTime = zeros(1,length(time_bins));
% rctC_trTime = zeros(1,length(time_bins));
% rct_sem_trTime = zeros(1,length(time_bins));
% rctC_sem_trTime = zeros(1,length(time_bins));
% rctE_trTime = zeros(1,length(time_bins));
% rctE_sem_trTime = zeros(1,length(time_bins));
nCyc_CRbase = nCyc_all;
nCyc_CRbase(eIx) = nCyc_all(eIx)-1;
CRbase_trL = [];
for itr = 1:length(nCyc_all)
   CRbase_trL = cat(2,CRbase_trL,double([cycTimeMs_all(itr):cycTimeMs_all(itr):cycTimeMs_all(itr)*nCyc_CRbase(itr)]));
end
failure_trL = bsxfun(@times,nCyc_all(eIx),cycTimeMs_all(eIx));
[h_time, bin_time_failures] = histc(failure_trL,trialTime_edges);
[h_time, bin_time_CRbase] = histc(CRbase_trL,trialTime_edges);

for ibin = 1:length(time_bins)
    ind_v = find(bin_time == time_bins(ibin));
    v = intersect(visInd_all,ind_v);    
    nOri_val_ind = cellfun(@intersect,ori_bin_indices,repmat({v},1,length(ori_bin_indices)),'unif',false);
    nOri_val = cell2mat(cellfun(@length, nOri_val_ind,'unif',false));
    
    c = intersect(cVisInd_all,ind_v);  
    nOri_inv_ind = cellfun(@intersect,oric_bin_indices,repmat({c},1,length(ori_bin_indices)),'unif',false);
    nOri_inv = cell2mat(cellfun(@length, nOri_inv_ind,'unif',false));
    
    if sum(lt(nOri_val,nOri_inv)) == 0
        multiplier = min(floor(nOri_val./nOri_inv));
        nOri_val_match = nOri_inv.*multiplier;
        nOri_inv_match = nOri_inv;
        ori_bin_match = cell(1,length(nOri_val));
        for i = 1:length(nOri_val)
            rand_ind = randperm(nOri_val(i),nOri_val_match(i));
            ori_bin_match{i} = nOri_val_ind{i}(rand_ind);
        end
        oric_bin_match = nOri_inv_ind;
    else
        error('cannot match data with existing code');
    end
    
    v = cell2mat(ori_bin_match);
    c = cell2mat(oric_bin_match);
    
    avg_trTime(ibin) = mean(trialTimeMs_all(ind_v));
%     rct_trTime(ibin) = mean(tarRct_all(intersect(find(sIx),v)));
%     rct_sem_trTime(ibin) = nanstd(tarRct_all(intersect(find(sIx),v)))/sqrt(length(intersect(find(sIx),v)));
%     rctC_trTime(ibin) = mean(cRct_all(intersect(find(faIx),c)));
%     rctC_sem_trTime(ibin) = nanstd(cRct_all(intersect(find(faIx),c)))/sqrt(length(intersect(find(faIx),c)));
%     rctE_trTime(ibin) = nanmean(eRct_all(intersect(find(eIx),v)));
%     rctE_sem_trTime(ibin) = nanstd(eRct_all(intersect(find(eIx),v)))/sqrt(length(intersect(find(eIx),v)));
    
    hits_trTime(ibin) = sum(sIx(v));
    misses_trTime(ibin) = sum(mIx(v));
    fa_trTime(ibin) = sum(faIx(c));
    cr_trTime(ibin) = sum(crIx(c));
    
    f = find(bin_time_failures == time_bins(ibin));
    crb = find(bin_time_CRbase == time_bins(ibin));
    earlies_trTime(ibin) = length(f);
    CRbase_trTime(ibin) = length(crb);
    
end

%min trial n indices
hrBinIx = find(hits_trTime+misses_trTime > minTrNall);
hBinIx = find(hits_trTime > minTrNall);
faBinIx = find(fa_trTime > minTrNall);
cBinIx = find(fa_trTime+cr_trTime > minTrNall);
eBinIx = find(earlies_trTime > minTrNall);

[HR_trTime ci_95_HR_trTime] = binofit(hits_trTime,hits_trTime+misses_trTime);
[FR_trTime ci_95_FR_trTime] = binofit(fa_trTime,fa_trTime+cr_trTime);
[ER_trTime ci_95_ER_trTime] = binofit(earlies_trTime,earlies_trTime+CRbase_trTime);

xTick_ind = find(~isnan(avg_trTime));

% workingMemoryFig = figure;
% subplot(2,1,1)
% errorbar(avg_trTime(hrBinIx),HR_trTime(hrBinIx),HR_trTime(hrBinIx)-ci_95_HR_trTime((hrBinIx),1)',ci_95_HR_trTime((hrBinIx),2)'-HR_trTime(hrBinIx),'go');
% hold on
% errorbar(avg_trTime(cBinIx),FR_trTime(cBinIx),FR_trTime(cBinIx)-ci_95_FR_trTime((cBinIx),1)',ci_95_FR_trTime((cBinIx),2)'-FR_trTime(cBinIx),'co');
% hold on
% errorbar(avg_trTime(eBinIx),ER_trTime(eBinIx),ER_trTime(eBinIx)-ci_95_ER_trTime((eBinIx),1)',ci_95_ER_trTime((eBinIx),2)'-ER_trTime(eBinIx),'mo');
% hold on
% xlim([0 trialTime_edges(end)+500])
% ylim([0 1])
% xlabel('trial length (s)')
% ylabel('HR')
% set(gca, 'XTick', avg_trTime(xTick_ind), ...
%        'XTickLabel', chop(avg_trTime(xTick_ind)/1000,2));
% axis square
% title('visual trials - HR')
% legend({'valid';'invalid';'earlies'},'Location','southeastoutside')

workingMemoryFig = figure;
subplot(2,1,1)
errorbar(avg_trTime,HR_trTime,HR_trTime-ci_95_HR_trTime(:,1)',ci_95_HR_trTime(:,2)'-HR_trTime,'go');
hold on
errorbar(avg_trTime,FR_trTime,FR_trTime-ci_95_FR_trTime(:,1)',ci_95_FR_trTime(:,2)'-FR_trTime,'co');
hold on
errorbar(avg_trTime,ER_trTime,ER_trTime-ci_95_ER_trTime(:,1)',ci_95_ER_trTime(:,2)'-ER_trTime,'mo');
hold on
xlim([0 trialTime_edges(end)+500])
ylim([0 1])
xlabel('trial length (s)')
ylabel('HR')
set(gca, 'XTick', avg_trTime(xTick_ind), ...
       'XTickLabel', chop(avg_trTime(xTick_ind)/1000,2));
axis square
title('visual trials - HR')
legend({'valid';'invalid';'earlies'},'Location','southeastoutside')

% figure(workingMemoryFig)
% subplot(2,1,2)
% errorbar(avg_trTime(hrBinIx) ,rct_trTime(hrBinIx) ,rct_sem_trTime(hrBinIx),'go')
% hold on
% errorbar(avg_trTime(faBinIx),rctC_trTime(faBinIx) ,rctC_sem_trTime(faBinIx),'co')
% hold on
% errorbar(avg_trTime(eBinIx),rctE_trTime(eBinIx),rctE_sem_trTime(eBinIx),'mo')
% hold on
% xlabel('trial length (ms)')
% ylabel('RT')
% xlim([0 trialTime_edges(end)+500])
% ylim([0 550])
% axis square
% title('visual trials - RT')
% legend({'valid';'invalid';'earlies'})

%% plot HR - vis trials, valid and invalid, for each cycle.

%% save 'uni bins' figs
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
set(0,'DefaultaxesFontSize', 16)

figure(uniHRFig);
print(fullfile(rc.fitOutputSummary, [date '_uniHR_fsavImg.pdf']),'-dpdf','-fillpage')
figure(uniRTFig);
print(fullfile(rc.fitOutputSummary, [date '_uniRT_fsavImg.pdf']),'-dpdf','-fillpage')
figure(workingMemoryFig);
print(fullfile(rc.fitOutputSummary, [date '_visbyTrialLength_fsavImg.pdf']),'-dpdf','-fillpage')
% figure(rtAll);
% print(fullfile(rc.fitOutputSummary, [date '_meanRT_fsavImg.pdf']),'-dpdf','-fillpage')
