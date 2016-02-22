clear all
close all
earlyLT50 = 1;
lapseLT10 = 1;
amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
av = behavParamsAV;
mice = unique(xd.Subject);
nMice = length(mice);

fn = fullfile(rc.fitOutputSummary, [date '_CatchSummary.mat']);
% fn = fullfile(rc.fitOutputSummary, ['17-Nov-2015_i613_i614_i616_CatchSummary.mat']);
load(fn)

AFig = figure;
VFig = figure;

for imouse = 1:nMice;
    mouse_name = mice(imouse);
    if earlyLT50 
        early_ind = find(mouse(imouse).early_mat<0.5);
    else
        early_ind = 1:size(mouse(imouse).early_mat,1);
    end
    if lapseLT10 
        lapse_ind = intersect(find(mouse(imouse).HR_ori_mat>0.9), find(mouse(imouse).HR_amp_mat>0.9));
    else
        lapse_ind = 1:size(mouse(imouse).HR_ori_mat,1);
    end
    
    use_ind = intersect(early_ind,lapse_ind);
    input = mouse(imouse).input(use_ind);
    input = concatenateStructures(input');
    
    missedIx = strcmp(input.trialOutcomeCell, 'ignore');
    successIx = strcmp(input.trialOutcomeCell, 'success');
    FAIx = strcmp(input.catchTrialOutcomeCell, 'FA');
    CRIx = strcmp(input.catchTrialOutcomeCell, 'CR');

    catchDirectionDeg = chop(double(celleqel2mat_padded(input.tCatchGratingDirectionDeg)),2);
    catchAmplitude = chop(double(celleqel2mat_padded(input.tSoundCatchAmplitude, NaN, 'double')),2);
    catchOris = unique(catchDirectionDeg);
    catchAmps = unique(catchAmplitude(~isnan(catchAmplitude)));
    targetDirectionDeg = chop(double(celleqel2mat_padded(input.tGratingDirectionDeg)),2);
    targetAmplitude = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude, NaN, 'double')),2);
    targetOris = unique(targetDirectionDeg);
    targetAmps = unique(targetAmplitude(~isnan(targetAmplitude)));
    
    targetReact = celleqel2mat_padded(input.reactTimeMs);
    catchReact = celleqel2mat_padded(input.leverUpTimeMs) - celleqel2mat_padded(input.tCatchTimeMs);
    
    % make cell array of target and catch values for success, miss, FA, and
    % CR trials
    tarRct{imouse} = targetReact;
    cRct{imouse} = catchReact;
    
    tarDeg{imouse} = targetDirectionDeg;
    tarAmp{imouse} = targetAmplitude;
    cDeg{imouse} = catchDirectionDeg;
    cAmp{imouse} = catchAmplitude;
    
    sIx_all{imouse} = successIx;
    mIx_all{imouse} = missedIx;
    faIx_all{imouse} = FAIx;
    crIx_all{imouse} = CRIx;
    

    [h_amp, bin_amp] = histc(targetAmplitude, amp_edges);
    [h_ampc, bin_ampc] = histc(catchAmplitude, amp_edges);
    [h_ori, bin_ori] = histc(targetDirectionDeg, ori_edges);
    [h_oric, bin_oric] = histc(catchDirectionDeg, ori_edges);

    amp_bins = unique(bin_amp);
    avg_amp = zeros(1,max(amp_bins,[],2));
    sem_amp = zeros(1,max(amp_bins,[],2));
    Hits_amp = zeros(1,max(amp_bins,[],2));
    Misses_amp = zeros(1,max(amp_bins,[],2));
    avg_ampc = zeros(1,max(amp_bins,[],2));
    sem_ampc = zeros(1,max(amp_bins,[],2));
    FAs_amp = zeros(1,max(amp_bins,[],2));
    CRs_amp = zeros(1,max(amp_bins,[],2));
    rct_amp = zeros(1,max(amp_bins,[],2));
    rctsem_amp = zeros(1,max(amp_bins,[],2));
    rct_ampc = zeros(1,max(amp_bins,[],2));
    rctsem_ampc = zeros(1,max(amp_bins,[],2));
    n_amp = zeros(2,max(amp_bins,[],2));
    for ibin = 1:max(amp_bins,[],2)
        ind = find(bin_amp == ibin);
        if length(ind)>5
            avg_amp(ibin) = mean(targetAmplitude(:,ind),2);
            sem_amp(ibin) = std(targetAmplitude,[],2)./sqrt(length(ind));
            Hits_amp(ibin) = sum(successIx(ind),2);
            Misses_amp(ibin) = sum(missedIx(ind),2);
            rct_amp(ibin) = mean(targetReact(successIx(ind)),2);
            rctsem_amp(ibin) = std(targetReact(successIx(ind)) ,[],2)./length(ind);
        end
        indc = find(bin_ampc == ibin);
        if length(indc)>5
            avg_ampc(ibin) = mean(catchAmplitude(:,indc),2);
            sem_ampc(ibin) = std(catchAmplitude(:,indc),[],2)./sqrt(length(indc));
            FAs_amp(ibin) = sum(FAIx(indc),2);
            CRs_amp(ibin) = sum(CRIx(indc),2);
            rct_ampc(ibin) = mean(catchReact(FAIx(indc)),2);
            rctsem_ampc(ibin) = std(catchReact(FAIx(indc)),[],2)./length(indc);
        end
    end
    n_amp = [Hits_amp+Misses_amp; FAs_amp+CRs_amp];
    [HR_amp, ci95_HR_amp] = binofit(Hits_amp, Misses_amp + Hits_amp);
    [FR_amp, ci95_FR_amp] = binofit(FAs_amp, CRs_amp + FAs_amp);

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
        if length(ind)>5
            avg_ori(ibin) = mean(targetDirectionDeg(:,ind),2);
            sem_ori(ibin) = std(targetDirectionDeg,[],2)./sqrt(length(ind));
            Hits_ori(ibin) = sum(successIx(ind),2);
            Misses_ori(ibin) = sum(missedIx(ind),2);
            rct_ori(ibin) = mean(targetReact(successIx(ind)),2)
            rctsem_ori(ibin) = std(targetReact(successIx(ind)),[],2)./length(ind);
        end
        indc = find(bin_oric == ibin);
        if length(indc)>5
            avg_oric(ibin) = mean(catchDirectionDeg(:,indc),2);
            sem_oric(ibin) = std(catchDirectionDeg(:,indc),[],2)./sqrt(length(indc));
            FAs_ori(ibin) = sum(FAIx(indc),2);
            CRs_ori(ibin) = sum(CRIx(indc),2);
            rct_oric(ibin) = mean(catchReact(FAIx(indc)),2);
            rctsem_oric(ibin) = std(catchReact(FAIx(indc)),[],2)./length(indc);
        end
    end
    [HR_ori, ci95_HR_ori] = binofit(Hits_ori, Misses_ori + Hits_ori);
    [FR_ori, ci95_FR_ori] = binofit(FAs_ori, CRs_ori + FAs_ori);
    n_ori = [Hits_ori+Misses_ori; FAs_ori+CRs_ori];
    figure; 
    errorbarxy(avg_amp, HR_amp, sem_amp, sem_amp, HR_amp - ci95_HR_amp(:,1)', ci95_HR_amp(:,2)' - HR_amp, {'ok', 'k', 'k'});
    hold on
    errorbarxy(avg_ampc, FR_amp, sem_ampc, sem_ampc, FR_amp - ci95_FR_amp(:,1)', ci95_FR_amp(:,2)' - FR_amp, {'oc', 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['i' num2str(mouse_name) '- Auditory and Auditory catch trials'])
    xlabel('Tone volume')
    ylabel('Hit rate')
    ylim([0 1])
    legend(['n = ' num2str(sum(n_amp(1,:),2))], ['n = ' num2str(sum(n_amp(2,:),2))]);
%     print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_AuditoryCatchSummary.pdf']),'-dpdf')

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
%     print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_VisualCatchSummary.pdf']),'-dpdf')
    
    figure(AFig);
    ind_catch = find(n_amp(2,:)>5);
    HitsA = zeros(1,length(ind_catch));
    MissesA = zeros(1,length(ind_catch));
    FAsA = zeros(1,length(ind_catch));
    CRsA = zeros(1,length(ind_catch));
    for ibin = 1:length(ind_catch)
       ampcs = unique(catchAmplitude(find(bin_ampc==(ind_catch(ibin)))));
       indT = find(ismember(targetAmplitude,ampcs));
       indC = find(ismember(catchAmplitude,ampcs));
       HitsA(ibin) = sum(successIx(indT),2);
       MissesA(ibin) = sum(missedIx(indT),2);
       FAsA(ibin) = sum(FAIx(indC),2);
       CRsA(ibin) = sum(CRIx(indC),2); 
    end
    try
    [HR_A, HR_ci95_A] = binofit(HitsA, HitsA+MissesA);
    [FR_A, FR_ci95_A] = binofit(FAsA, FAsA+CRsA);
    errorbarxy(HR_A, FR_A, HR_A-HR_ci95_A(:,1)', HR_ci95_A(:,2)'-HR_A, FR_A-FR_ci95_A(:,1)', FR_ci95_A(:,2)'-FR_A, {['o' av(imouse).col_str], av(imouse).col_str, av(imouse).col_str})
    hold on
    x = 0:.1:1;
    y = x;
    plot(x,y,'--k')
    hold on
    xlim([0 1])
    ylim([0 1])
    xlabel('Hit rate- valid cue')
    ylabel('Hit rate- invalid cue')
    title('Auditory trials')
    catch
        disp('no aud trials to plot')
    end
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
    HitsA = zeros(1,length(ind_catch));
    MissesA = zeros(1,length(ind_catch));
    FAsA = zeros(1,length(ind_catch));
    CRsA = zeros(1,length(ind_catch));
    if ~isempty(ind_catch)
    for ibin = 1:length(ind_catch)
       ampcs = unique(catchDirectionDeg(find(bin_oric==(ind_catch(ibin)))));
       indT = find(ismember(targetDirectionDeg,ampcs));
       indC = find(ismember(catchDirectionDeg,ampcs));
       HitsV(ibin) = sum(successIx(indT),2);
       MissesV(ibin) = sum(missedIx(indT),2);
       FAsV(ibin) = sum(FAIx(indC),2);
       CRsV(ibin) = sum(CRIx(indC),2); 
    end
    [HR_V, HR_ci95_V] = binofit(HitsV, HitsV+MissesV);
    [FR_V, FR_ci95_V] = binofit(FAsV, FAsV+CRsV);
    errorbarxy(HR_V, FR_V, HR_V-HR_ci95_V(:,1)', HR_ci95_V(:,2)'-HR_V, FR_V-FR_ci95_V(:,1)', FR_ci95_V(:,2)'-FR_V, {['o' av(imouse).col_str], av(imouse).col_str, av(imouse).col_str})
    hold on
    xlim([0 1])
    ylim([0 1])
    plot(x,y,'--k')
    hold on
    xlabel('Hit rate- valid cue')
    ylabel('Hit rate- invalid cue')
    title('Visual trials')
    end

end

figure(VFig);
print(fullfile(rc.fitOutputSummary, [date '_i613_i614_i616_VisualCatchSummary.pdf']),'-dpdf')
figure(AFig);
print(fullfile(rc.fitOutputSummary, [date '_i613_i614_i616_AuditoryCatchSummary.pdf']),'-dpdf')


%% plot HR summary across mice

% create universal ori and amp bins
tarDeg_all = cat(2,cell2mat(tarDeg));
cDeg_all = cat(2,cell2mat(cDeg));
tarAmp_all = cat(2,cell2mat(tarAmp));
cAmp_all = cat(2,cell2mat(cAmp));
sIx = cat(2,cell2mat(sIx_all));
mIx = cat(2,cell2mat(mIx_all));
faIx = cat(2,cell2mat(faIx_all));
crIx = cat(2,cell2mat(crIx_all));
tarRct_all = cat(2,cell2mat(tarRct));
cRct_all = cat(2,cell2mat(cRct));

reducedOriC_edges = [1 24 100];
catch_ori_edges = unique(cDeg_all(~isnan(cDeg_all)));
catch_amp_edges = unique(cAmp_all(~isnan(cAmp_all)));

[h_ori, bin_ori_all] = histc(tarDeg_all, ori_edges);
[h_oric, bin_oric_all] = histc(cDeg_all, ori_edges);
[h_amp, bin_amp_all] = histc(tarAmp_all, amp_edges);
[h_ampc, bin_ampc_all] = histc(cAmp_all, amp_edges);

ori_bins = unique(bin_ori_all);
amp_bins = unique(bin_amp_all);
avg_ori_all = zeros(1,max(ori_bins,[],2));
sem_ori_all = zeros(1,max(ori_bins,[],2));    
avg_amp_all = zeros(1,max(amp_bins,[],2));
sem_amp_all = zeros(1,max(amp_bins,[],2));
avg_oric_all = zeros(1,max(ori_bins,[],2));
sem_oric_all = zeros(1,max(ori_bins,[],2));    
avg_ampc_all = zeros(1,max(amp_bins,[],2));
sem_ampc_all = zeros(1,max(amp_bins,[],2));
rct_ori_all = NaN(1,max(ori_bins,[],2));
rctsem_ori_all = NaN(1,max(ori_bins,[],2));
rct_oric_all = NaN(1,max(ori_bins,[],2));
rctsem_oric_all = NaN(1,max(ori_bins,[],2));
rct_amp_all = NaN(1,max(amp_bins,[],2));
rctsem_amp_all = NaN(1,max(amp_bins,[],2));
rct_ampc_all = NaN(1,max(amp_bins,[],2));
rctsem_ampc_all = NaN(1,max(amp_bins,[],2));
for ibin = 1:max(ori_bins,[],2)
    ind = find(bin_ori_all == ibin);
    avg_ori_all(ibin) = mean(tarDeg_all(ind),2);
    sem_ori_all(ibin) = std(tarDeg_all(ind),[],2)./length(ind);
    rct_ori_all(ibin) = mean(tarRct_all(intersect(find(sIx),ind)),2);
    rctsem_ori_all(ibin) = std(tarRct_all(intersect(find(sIx),ind)),[],2)./length(ind);
    indc = find(bin_oric_all == ibin);
    avg_oric_all(ibin) = mean(cDeg_all(indc),2);
    sem_oric_all(ibin) = std(cDeg_all(indc),[],2)./length(indc);
    rct_oric_all(ibin) = mean(cRct_all(intersect(find(faIx),indc)),2);
    rctsem_oric_all(ibin) = std(cRct_all(intersect(find(faIx),indc)),[],2)./length(indc);
end

for ibin = 1:max(amp_bins,[],2);
    ind = find(bin_amp_all == ibin);
    avg_amp_all(ibin) = mean(tarAmp_all(ind),2);
    sem_amp_all(ibin) = std(tarAmp_all(ind),[],2)./length(ind);
    rct_amp_all = mean(tarRct_all(intersect(find(sIx),ind)),2);
    rctsem_amp_all = std(tarRct_all(intersect(find(sIx),ind)),[],2)./length(ind);
    indc = find(bin_ampc_all == ibin);
    avg_ampc_all(ibin) = mean(cAmp_all(indc),2);
    sem_ampc_all(ibin) = std(cAmp_all(indc),[],2)./length(indc);
    rct_ampc_all = mean(cRct_all(intersect(find(faIx),indc)),2);
    rctsem_ampc_all =  std(cRct_all(intersect(find(faIx),indc)),[],2)./length(indc);
end


hitsAll_V = zeros(1,max(ori_bins,[],2));
missAll_V = zeros(1,max(ori_bins,[],2));
faAll_V = zeros(1,max(ori_bins,[],2));
crAll_V = zeros(1,max(ori_bins,[],2));
n_ori = zeros(1,max(ori_bins,[],2));
for ibin = 1:max(ori_bins,[],2)
    ind = find(bin_ori_all == ibin);
    if sum(sIx(ind)) + sum(mIx(ind)) > 10
    hitsAll_V(ibin) = sum(sIx(ind));
    missAll_V(ibin)= sum(mIx(ind));
    end
    indc = find(bin_oric_all == ibin);
    if sum(faIx(indc)) + sum(crIx(indc)) > 10
    faAll_V(ibin) = sum(faIx(indc));
    crAll_V(ibin) = sum(crIx(indc));
    end
    n_ori(ibin) = hitsAll_V(ibin)+missAll_V(ibin);
    n_oric(ibin) = faAll_V(ibin)+crAll_V(ibin);
end

[HR_ori_all, ci95_HR_ori_all] = binofit(hitsAll_V, missAll_V + hitsAll_V);
[FR_ori_all, ci95_FR_ori_all] = binofit(faAll_V, faAll_V + crAll_V);

hitsAll_A = zeros(1,max(amp_bins,[],2));
missAll_A= zeros(1,max(amp_bins,[],2));
faAll_A = zeros(1,max(amp_bins,[],2));
crAll_A= zeros(1,max(amp_bins,[],2));
for ibin = 1:max(amp_bins,[],2)
    ind = find(bin_amp_all == ibin);
    if  sum(sIx(ind)) + sum(mIx(ind)) >10
    hitsAll_A(ibin) = sum(sIx(ind));
    missAll_A(ibin)= sum(mIx(ind));
    end
    indc = find(bin_ampc_all == ibin);
    if  sum(faIx(indc)) + sum(crIx(indc)) >10
    faAll_A(ibin) = sum(faIx(indc));
    crAll_A(ibin) = sum(crIx(indc));
    end
    n_amp(ibin) = hitsAll_A(ibin)+missAll_A(ibin);
    n_ampc(ibin) = faAll_A(ibin)+crAll_A(ibin);
end

[HR_amp_all, ci95_HR_amp_all] = binofit(hitsAll_A, missAll_A + hitsAll_A);
[FR_amp_all, ci95_FR_amp_all] = binofit(faAll_A, faAll_A + crAll_A);


nTrials = cell2mat(cellfun(@length,tarDeg,'Unif',false));
cumTrials = cumsum(nTrials);

%%
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
    
    Hits_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
    Misses_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
    FAs_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
    CRs_amp_all{imouse} = zeros(1,max(amp_bins,[],2));
    
    n_ori_all{imouse} = zeros(2,max(amp_bins,[],2));
    n_amp_all{imouse} = zeros(2,max(amp_bins,[],2));
    
    rct_ori_each{imouse} = NaN(1,max(ori_bins,[],2));
    rctsem_ori_each{imouse} = NaN(1,max(ori_bins,[],2));
    rct_oric_each{imouse} = NaN(1,max(ori_bins,[],2));
    rctsem_oric_each{imouse} = NaN(1,max(ori_bins,[],2));
    rct_amp_each{imouse} = NaN(1,max(amp_bins,[],2));
    rctsem_amp_each{imouse} = NaN(1,max(amp_bins,[],2));
    rct_ampc_each{imouse} = NaN(1,max(amp_bins,[],2));
    rctsem_ampc_each{imouse} = NaN(1,max(amp_bins,[],2));
    
    for ibin = 1:max(ori_bins,[],2)
        ind = find(bin_ori_all(trialInd) == ibin);
        if length(ind)>5
            Hits_ori_all{imouse}(ibin) = sum(sIx_all{imouse}(ind),2);
            Misses_ori_all{imouse}(ibin) = sum(mIx_all{imouse}(ind),2);
            rct_ori_each{imouse}(ibin) = mean(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind)),2);
            rctsem_ori_each{imouse}(ibin) = std(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind)),[],2)./length(ind);
        end
        indc = find(bin_oric_all(trialInd) == ibin);
        if length(indc)>5
            FAs_ori_all{imouse}(ibin) = sum(faIx_all{imouse}(indc),2);
            CRs_ori_all{imouse}(ibin) = sum(crIx_all{imouse}(indc),2);
            rct_oric_each{imouse}(ibin) = mean(cRct{imouse}(intersect(find(faIx_all{imouse}),indc)),2);
            rctsem_oric_each{imouse}(ibin) = std(cRct{imouse}(intersect(find(faIx_all{imouse}),indc)),[],2)./length(indc);
        end
    end
    [HR_ori_each{imouse}, ci95_HR_ori_each{imouse}] = binofit(Hits_ori_all{imouse}, Misses_ori_all{imouse} + Hits_ori_all{imouse});
    [FR_ori_each{imouse}, ci95_FR_ori_each{imouse}] = binofit(FAs_ori_all{imouse}, CRs_ori_all{imouse} + FAs_ori_all{imouse});
    n_ori_all{imouse} = [Hits_ori_all{imouse}+Misses_ori_all{imouse}; FAs_ori_all{imouse}+CRs_ori_all{imouse}];

    for ibin = 1:max(amp_bins,[],2)
        ind = find(bin_amp_all(trialInd) == ibin);
        if length(ind)>5
            Hits_amp_all{imouse}(ibin) = sum(sIx_all{imouse}(ind),2);
            Misses_amp_all{imouse}(ibin) = sum(mIx_all{imouse}(ind),2);
            rct_amp_each{imouse}(ibin) = mean(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind)),2);
            rctsem_amp_each{imouse}(ibin) = std(tarRct{imouse}(intersect(find(sIx_all{imouse}),ind)),[],2)./length(ind);
        end
        indc = find(bin_ampc_all(trialInd) == ibin);
        if length(indc)>5
            FAs_amp_all{imouse}(ibin) = sum(faIx_all{imouse}(indc),2);
            CRs_amp_all{imouse}(ibin) = sum(crIx_all{imouse}(indc),2);
            rct_ampc_each{imouse}(ibin) = mean(cRct{imouse}(intersect(find(faIx_all{imouse}),indc)),2);
            rctsem_ampc_each{imouse}(ibin) = std(cRct{imouse}(intersect(find(faIx_all{imouse}),indc)),[],2)./length(indc);
        end
    end
    [HR_amp_each{imouse}, ci95_HR_amp_each{imouse}] = binofit(Hits_amp_all{imouse}, Misses_amp_all{imouse} + Hits_amp_all{imouse});
    [FR_amp_each{imouse}, ci95_FR_amp_each{imouse}] = binofit(FAs_amp_all{imouse}, CRs_amp_all{imouse} + FAs_amp_all{imouse});
    n_amp_all{imouse} = [Hits_amp_all{imouse}+Misses_amp_all{imouse}; FAs_amp_all{imouse}+CRs_amp_all{imouse}];
end

%% plot each mouse curve and avg across mice with universal bins
uniHRFig_V = figure;
uniHRFig_A = figure;
uniRTFig_V = figure;
uniRTFig_A = figure;

%visual trials
for imouse = 1:nMice
    figure(uniHRFig_V)
    hold on
    colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
    
    errorbarxy(avg_ori_all, HR_ori_each{imouse}, sem_ori_all, sem_ori_all, HR_ori_each{imouse} - ci95_HR_ori_each{imouse}(:,1)', ci95_HR_ori_each{imouse}(:,2)' - HR_ori_each{imouse}, {['o' 'k'], 'k', 'k'});
    hold on
    errorbarxy(avg_oric_all, FR_ori_each{imouse}, sem_oric_all, sem_oric_all, FR_ori_each{imouse} - ci95_FR_ori_each{imouse}(:,1)', ci95_FR_ori_each{imouse}(:,2)' - FR_ori_each{imouse}, {['o' 'c'], 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['visual trials'])
    xlabel('Orientation change (deg)')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([0 100])
%     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
    
end
figure(uniHRFig_V)
hold on
errorbarxy(avg_ori_all, HR_ori_all, sem_ori_all, sem_ori_all, HR_ori_all - ci95_HR_ori_all(:,1)', ci95_HR_ori_all(:,2)' - HR_ori_all, {['o' 'k'], 'k', 'k'});
hold on
errorbarxy(avg_oric_all, FR_ori_all, sem_oric_all, sem_oric_all, FR_ori_all - ci95_FR_ori_all(:,1)', ci95_FR_ori_all(:,2)' - FR_ori_all, {['o' 'c'], 'c', 'c'});
hold on
scatter(avg_ori_all,HR_ori_all,'ko','filled');
hold on
scatter(avg_oric_all,FR_ori_all,'co','filled')

%visual trials - RT
for imouse = 1:nMice
    figure(uniRTFig_V)
    hold on
    colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
    
    errorbarxy(avg_ori_all, rct_ori_each{imouse}, sem_ori_all, sem_ori_all, rct_ori_each{imouse} - rctsem_ori_each{imouse}, rctsem_ori_each{imouse} - rct_ori_each{imouse}, {['o' 'k'], 'k', 'k'});
    hold on
    errorbarxy(avg_oric_all, rct_oric_each{imouse}, sem_oric_all, sem_oric_all, rct_oric_each{imouse} - rctsem_ori_each{imouse}, rctsem_oric_each{imouse} - rct_oric_each{imouse}, {['o' 'c'], 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['visual trials'])
    xlabel('Orientation change (deg)')
    ylabel('RT')
%     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
    
end
figure(uniRTFig_V)
hold on
errorbarxy(avg_ori_all, rct_ori_all, sem_ori_all, sem_ori_all, rct_ori_all - rctsem_ori_all, rctsem_ori_all - rct_ori_all, {['o' 'k'], 'k', 'k'});
hold on
errorbarxy(avg_oric_all, rct_oric_all, sem_oric_all, sem_oric_all, rct_oric_all - rctsem_oric_all, rctsem_oric_all - rct_oric_all, {['o' 'c'], 'c', 'c'});
hold on
scatter(avg_ori_all,rct_ori_all,'ko','filled');
hold on
scatter(avg_oric_all,rct_oric_all,'co','filled');
hold on

%auditory trials
for imouse = 1:nMice
    figure(uniHRFig_A)
    hold on
    colstr = av(find(cell2mat({av.mouse}) == mice(imouse))).col_str;
    
    errorbarxy(avg_amp_all, HR_amp_each{imouse}, sem_amp_all, sem_amp_all, HR_amp_each{imouse} - ci95_HR_amp_each{imouse}(:,1)', ci95_HR_amp_each{imouse}(:,2)' - HR_amp_each{imouse}, {['o' 'k'], 'k', 'k'});
    hold on
    errorbarxy(avg_ampc_all, FR_amp_each{imouse}, sem_ampc_all, sem_ampc_all, FR_amp_each{imouse} - ci95_FR_amp_each{imouse}(:,1)', ci95_FR_amp_each{imouse}(:,2)' - FR_amp_each{imouse}, {['o' 'c'], 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['auditory trials'])
    xlabel('Target Sound Amplitude')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([0 1])
%     legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
    
end
figure(uniHRFig_A)
hold on
errorbarxy(avg_amp_all, HR_amp_all, sem_amp_all, sem_amp_all, HR_amp_all - ci95_HR_amp_all(:,1)', ci95_HR_amp_all(:,2)' - HR_amp_all, {['o' 'k'], 'k', 'k'});
hold on
errorbarxy(avg_ampc_all, FR_amp_all, sem_ampc_all, sem_ampc_all, FR_amp_all - ci95_FR_amp_all(:,1)', ci95_FR_amp_all(:,2)' - FR_amp_all, {['o' 'c'], 'c', 'c'});
hold on
scatter(avg_amp_all,HR_amp_all,'ko','filled');
hold on
scatter(avg_ampc_all,FR_amp_all,'co','filled')
    