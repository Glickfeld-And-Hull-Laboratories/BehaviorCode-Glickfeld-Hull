earlyLT50 = 1;
lapseLT10 = 1;
amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
miceAnalyzed = unique(xd.Subject);
av = behavParamsAV;
fn = fullfile(rc.fitOutputSummary, [date '_i613_i614_CatchSummary.mat']);
% fn = fullfile(rc.fitOutputSummary, ['11-Oct-2015_i613_i614_CatchSummary.mat']);
load(fn)

close all

AFig = figure;
VFig = figure;
%need to change bc need to turn off thresholds for 2x and 4x
maxEarlyRate = 1.0;
maxLapseRate = 0.0;

%for graphs, need to manually change titles

for imouse = 1:length(miceAnalyzed);
    mouse_name = av(imouse).mouse; 
    if earlyLT50 
        early_ind = find(mouse(imouse).early_mat<maxEarlyRate);
    else
        early_ind = 1:size(early_mat,1);
    end
    if lapseLT10 
        lapse_ind = intersect(find(mouse(imouse).HR_ori_mat>maxLapseRate), find(mouse(imouse).HR_amp_mat>maxLapseRate));
    else
        lapse_ind = 1:size(HR_ori_mat,1);
    end
    
    use_ind = intersect(early_ind,lapse_ind);
    input = mouse(imouse).input(use_ind);
    input = concatenateStructures(input');
    
     missedIx = strcmp(input.trialOutcomeCell, 'ignore');
     successIx = strcmp(input.trialOutcomeCell, 'success');
     FAIx = strcmp(input.catchTrialOutcomeCell, 'FA');
     
     CRIx = strcmp(input.catchTrialOutcomeCell, 'CR');
    
    
    catchDirectionDeg = chop(double(cell2mat_padded(input.tCatchGratingDirectionDeg)),2)';
   
       catchAmplitude = chop(double(celleqel2mat_padded(input.tSoundCatchAmplitude, NaN, 'double')),2);
% catchAmplitude = celleqel2mat_padded(input.tSoundCatchAmplitude);
    oriindx = find(catchDirectionDeg>0);
    ampindx = find(catchAmplitude>0);
    catchOris = unique(catchDirectionDeg);
    catchAmps = unique(catchAmplitude(~isnan(catchAmplitude)));
    targetDirectionDeg = chop(double(cell2mat_padded(input.tGratingDirectionDeg)),2)';
    targetAmplitude = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude, NaN, 'double')),2);
    targetOris = unique(targetDirectionDeg);
    targetAmps = unique(targetAmplitude(~isnan(targetAmplitude)));
    
    audtrials = find(cell2mat_padded(input.tGratingDirectionDeg) == 0);
    vistrials = find(cell2mat_padded(input.tGratingDirectionDeg) >0);
    
    successes = find(successIx == 1);
    FAs = find(FAIx ==1);
    successaud = intersect(successes, audtrials);
    successvis = intersect(successes, vistrials);
    
    FAaud = intersect(FAs, vistrials);
    FAvis = intersect(FAs, audtrials);
    
    
    
    reactV = celleqel2mat_padded(input.reactTimeMs);
    reactVsuccess = reactV(successIx);
    reactVsuccessa = reactV(successaud);
    reactVsuccessv= reactV(successvis);
%    reactC =  celleqel2mat_padded(input.leverUpTimeMs) -( celleqel2mat_padded(input.leverDownTimeMs) +((unique(input.stimOnTimeMs)+unique(input.stimOffTimeMs)).*(celleqel2mat_padded(input.catchCyclesOn)-1)));
    reactC = celleqel2mat_padded(input.leverUpTimeMs) - celleqel2mat_padded(input.tCatchTimeMs);
    reactCsuccess = reactC(successIx);
    reactCsuccessa = reactC(FAaud);
    reactCsuccessv = reactC(FAvis);
    targetsuccess = targetDirectionDeg(successvis);
    catchFA = catchDirectionDeg(FAvis);
    atargetsuccess = targetAmplitude(successaud);
    acatchFA = catchAmplitude(FAaud);
    
    
    
    
    [h_amp, bin_amp] = histc(atargetsuccess, amp_edges);%count how many things are in the bins
    [h_ampc, bin_ampc] = histc(acatchFA, amp_edges);
    [h_ori, bin_ori] = histc(targetsuccess, ori_edges);
    [h_oric, bin_oric] = histc(catchFA, ori_edges);

    amp_bins = unique(bin_amp);
    avg_amp = zeros(1,max(amp_bins,[],2));
    sem_amp = zeros(1,max(amp_bins,[],2));
    Hits_amp = zeros(1,max(amp_bins,[],2));
    Misses_amp = zeros(1,max(amp_bins,[],2));
    avg_ampc = zeros(1,max(amp_bins,[],2));
    sem_ampc = zeros(1,max(amp_bins,[],2));
    FAs_amp = zeros(1,max(amp_bins,[],2));
    CRs_amp = zeros(1,max(amp_bins,[],2));
    n_amp = zeros(2,max(amp_bins,[],2));
    success_amp = zeros(1,max(amp_bins,[],2));
    success_ampc = zeros(1,max(amp_bins,[],2));
    reacttime_amp = zeros(1,max(amp_bins,[],2));
     reacttime_ampc = zeros(1,max(amp_bins,[],2));
     reacttime_sem_amp = zeros(1,max(amp_bins,[],2));
     reacttime_sem_ampc = zeros(1,max(amp_bins,[],2));
    for ibin = 1:max(amp_bins,[],2)
        ind = find(bin_amp == ibin);
        if length(ind)>5%whyfive?
            avg_amp(ibin) = mean(atargetsuccess(:,ind),2);
            sem_amp(ibin) = std(atargetsuccess,[],2)./sqrt(length(ind));
            Hits_amp(ibin) = sum(successIx(ind),2);
            Misses_amp(ibin) = sum(missedIx(ind),2);
            success_amp(ibin) = length(ind);
            reacttime_amp(ibin) = mean(reactVsuccessa(ind),2);
            reacttime_sem_amp(ibin) = std(reactVsuccessa(ind),[],2)/sqrt(length(reactVsuccessa(ind)));
        end
        indc = find(bin_ampc == ibin);
        if length(indc)>5
            avg_ampc(ibin) = mean(acatchFA(:,indc),2);
            sem_ampc(ibin) = std(acatchFA,[],2)./sqrt(length(indc));
            FAs_amp(ibin) = sum(FAIx(indc),2);
            CRs_amp(ibin) = sum(CRIx(indc),2);
            success_ampc(ibin) = length(indc);
            reacttime_ampc(ibin) = mean(reactCsuccessa(indc),2);
            reacttime_sem_ampc(ibin) = std(reactCsuccessa(indc),[],2)/sqrt(length(reactCsuccessa(indc)));
        end
    end
    

 
%      n_amp = [Hits_amp+Misses_amp; FAs_amp+CRs_amp];
%     [HR_amp, ci95_HR_amp] = binofit(Hits_amp, Misses_amp + Hits_amp);
%     [FR_amp, ci95_FR_amp] = binofit(FAs_amp, CRs_amp + FAs_amp);

    ori_bins = unique(bin_ori);
    avg_ori = zeros(1,max(ori_bins,[],2));
    sem_ori = zeros(1,max(ori_bins,[],2));
    Hits_ori = zeros(1,max(ori_bins,[],2));
    Misses_ori = zeros(1,max(ori_bins,[],2));
    avg_oric = zeros(1,max(ori_bins,[],2));
    sem_oric = zeros(1,max(ori_bins,[],2));
    FAs_ori = zeros(1,max(ori_bins,[],2));
    CRs_ori = zeros(1,max(ori_bins,[],2));
    success_ori = zeros(1,max(ori_bins,[],2));
    success_oric = zeros(1,max(ori_bins,[],2));
    n_ori = zeros(2,max(ori_bins,[],2));
     reacttime_ori = zeros(1,max(ori_bins,[],2));
     reacttime_oric = zeros(1,max(ori_bins,[],2));
     reacttime_sem_ori = zeros(1,max(ori_bins,[],2));
     reacttime_sem_oric = zeros(1,max(ori_bins,[],2));
    for ibin = 1:max(ori_bins,[],2)
        ind = find(bin_ori == ibin);
         if length(ind)>5
            avg_ori(ibin) = mean(targetsuccess(:,ind),2);
            sem_ori(ibin) = std(targetsuccess,[],2)./sqrt(length(ind));
            Hits_ori(ibin) = sum(successIx(ind),2);
            Misses_ori(ibin) = sum(missedIx(ind),2);
            success_ori(ibin) = length(ind);
            reacttime_ori(ibin) = mean(reactVsuccessv(ind),2);
            reacttime_sem_ori(ibin) = std(reactVsuccessv(ind),[],2)/sqrt(length(reactVsuccessv(ind)));
        end
        indc = find(bin_oric == ibin);
%          if length(indc)<5
            avg_oric(ibin) = mean(catchFA(:,indc),2);
            sem_oric(ibin) = std(catchFA,[],2)./sqrt(length(indc));
            FAs_ori(ibin) = sum(FAIx(indc),2);
            CRs_ori(ibin) = sum(CRIx(indc),2);
            success_oric(ibin) = length(indc);
            reacttime_oric(ibin) = mean(reactCsuccessv(indc),2);
            reacttime_sem_oric(ibin) = std(reactCsuccessv(indc),[],2)/sqrt(length(reactCsuccessv(indc)));
    
    end
        

 
%      n_amp = [Hits_amp+Misses_amp; FAs_amp+CRs_amp];
%     [HR_ori, ci95_HR_ori] = binofit(Hits_ori, Misses_ori + Hits_ori);
%     [FR_ori, ci95_FR_ori] = binofit(FAs_ori, CRs_ori + FAs_ori);
%      n_ori = [Hits_ori+Misses_ori; FAs_ori+CRs_ori];
        n_amp = [success_amp; success_ampc];
        
        n_ori = [success_ori;success_oric];




       figure; 
       errorbarxy(avg_amp, reacttime_amp, sem_amp, reacttime_sem_amp,{'ok', 'k', 'k'});
    hold on
    errorbarxy(avg_ampc, reacttime_ampc, sem_ampc, reacttime_sem_ampc, {'oc', 'c', 'c'});
    
%     errorbarxy(avg_amp, HR_amp, sem_amp, sem_amp, HR_amp - ci95_HR_amp(:,1)', ci95_HR_amp(:,2)' - HR_amp, {'ok', 'k', 'k'});
%     hold on
%     errorbarxy(avg_ampc, FR_amp, sem_ampc, sem_ampc, FR_amp - ci95_FR_amp(:,1)', ci95_FR_amp(:,2)' - FR_amp, {'oc', 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['i' num2str(mouse_name) '- Auditory and Auditory catch trials 2x'])
    xlabel('Tone volume')
    ylabel('React Time')
%     ylim([min(reacttime_amp) max(reacttime_ampc)+ 200])
%     xlim([-200 max(avg_ampc)+ 100])
    legend(['n = ' num2str(sum(n_amp(1,:),2))], ['n = ' num2str(sum(n_amp(2,:),2))]);
    print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_AuditoryCatchSummary.pdf']),'-dpdf')

    figure; 


    errorbarxy(avg_ori, reacttime_ori, sem_ori, reacttime_sem_ori,{'ok', 'k', 'k'});
    hold on
    errorbarxy(avg_oric, reacttime_oric, sem_oric, reacttime_sem_oric, {'oc', 'c', 'c'});
      
      
%     errorbarxy(avg_ori, HR_ori, sem_ori, sem_ori, HR_ori - ci95_HR_ori(:,1)', ci95_HR_ori(:,2)' - HR_ori, {'ok', 'k', 'k'});
%     hold on
%     errorbarxy(avg_oric, FR_ori, sem_oric, sem_oric, FR_ori - ci95_FR_ori(:,1)', ci95_FR_ori(:,2)' - FR_ori, {'oc', 'c', 'c'});
    set(gca, 'xscale', 'log');
    title(['i' num2str(mouse_name) '- Visual and visual catch trials 2x'])
    xlabel('Orientation change (deg)')
    ylabel('React Time')
     ylim([-200 max(reacttime_oric) + 500])
     xlim([0 100])
    legend(['n = ' num2str(sum(n_ori(1,:),2))], ['n = ' num2str(sum(n_ori(2,:),2))]);
    print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_VisualCatchSummary.pdf']),'-dpdf')
    
    figure(AFig);
 
    
    ind_catch = find(n_amp(2,:)>5);
    HitsA = zeros(1,length(ind_catch));
    MissesA = zeros(1,length(ind_catch));
    FAsA = zeros(1,length(ind_catch));
    CRsA = zeros(1,length(ind_catch));
    reactamp = zeros(1, length(ind_catch));
    reactampc = zeros(1,length(ind_catch));
    react_sem_amp= zeros(1,length(ind_catch));
    react_sem_ampc= zeros(1,length(ind_catch));
    for ibin = 1:length(ind_catch)
       ampcs = unique(catchAmplitude(find(bin_ampc==(ind_catch(ibin)))));
       indT = find(ismember(targetsuccess,ampcs));
       indC = find(ismember(acatchFA,ampcs));
       HitsA(ibin) = sum(successIx(indT),2);
       MissesA(ibin) = sum(missedIx(indT),2);
       FAsA(ibin) = sum(FAIx(indC),2);
       CRsA(ibin) = sum(CRIx(indC),2);
       reactamp(ibin) = mean(reactVsuccessa(indT),2);
       reactampc(ibin) =  mean(reactCsuccessa(indC),2);
       react_sem_amp(ibin) = std(reactVsuccessa(indT),[],2)/sqrt(length(reactVsuccessa(indT)));
       react_sem_ampc(ibin) = std(reactCsuccessa(indC),[],2)/sqrt(length(reactCsuccessa(indC)));
        
    end
%     [HR_A, HR_ci95_A] = binofit(HitsA, HitsA+MissesA);
%     [FR_A, FR_ci95_A] = binofit(FAsA, FAsA+CRsA);

    
    errorbarxy(reactamp, reactampc, react_sem_amp, react_sem_ampc); 
%     errorbarxy(HR_A, FR_A, HR_A-HR_ci95_A(:,1)', HR_ci95_A(:,2)'-HR_A, FR_A-FR_ci95_A(:,1)', FR_ci95_A(:,2)'-FR_A, {['o' av(imouse).col_str], av(imouse).col_str, av(imouse).col_str})
    hold on
%     x = 0:.1:1;
%     y = x;
%     plot(x,y,'--k')
    hold on
%     xlim([0 1])
%     ylim([0 1])
    xlabel('React time- valid cue')
    ylabel('React time- invalid cue')
    title('Auditory trials-2x')
    
    figure(VFig);
    ind_catch = find(n_ori(2,:)>1);
    HitsA = zeros(1,length(ind_catch));
    MissesA = zeros(1,length(ind_catch));
    FAsA = zeros(1,length(ind_catch));
    CRsA = zeros(1,length(ind_catch));
    reactori = zeros(1, length(ind_catch));
    reactoric = zeros(1,length(ind_catch));
    react_sem_ori= zeros(1,length(ind_catch));
    react_sem_oric= zeros(1,length(ind_catch));
    
    for ibin = 1:length(ind_catch)
       ampcs = unique(catchFA(find(bin_oric==(ind_catch(ibin)))));
       indT = find(ismember(targetsuccess,ampcs));
       indC = find(ismember(catchFA,ampcs));
       HitsV(ibin) = sum(successIx(indT),2);
       MissesV(ibin) = sum(missedIx(indT),2);
       FAsV(ibin) = sum(FAIx(indC),2);
       CRsV(ibin) = sum(CRIx(indC),2);
       reactori(ibin) = mean(reactVsuccessv(indT),2);
       reactoric(ibin) =  mean(reactCsuccessv(indC),2);
       react_sem_ori(ibin) = std(reactVsuccessv(indT),[],2)/sqrt(length(reactVsuccessv(indT)));
       react_sem_oric(ibin) = std(reactCsuccessv(indC),[],2)/sqrt(length(reactCsuccessv(indC)));
        
       
    end
%     [HR_V, HR_ci95_V] = binofit(HitsV, HitsV+MissesV);
%     [FR_V, FR_ci95_V] = binofit(FAsV, FAsV+CRsV);
   
    
    errorbarxy(reactori, reactoric, react_sem_ori, react_sem_oric, {['o' av(imouse).col_str], av(imouse).col_str, av(imouse).col_str}); 
%     errorbarxy(HR_V, FR_V, HR_V-HR_ci95_V(:,1)', HR_ci95_V(:,2)'-HR_V, FR_V-FR_ci95_V(:,1)', FR_ci95_V(:,2)'-FR_V, {['o' av(imouse).col_str], av(imouse).col_str, av(imouse).col_str})
    hold on
%     xlim([0,1000])
%     ylim([0,1000])
%     xlim([0 300])
%     ylim([0 300])
%      xlim([0 2000])
%      ylim([0 2000])
    hold on
    x = 0:.1: max(reactoric)+200;
    y = x;
    plot(x,y,'--k')
    hold on
   
    xlabel('React time- valid cue')
    ylabel('React time- invalid cue')
    title('Visual trials -2x')
     end
figure(VFig);
% print(fullfile(rc.fitOutputSummary, [date '_i613_i614_VisualCatchSummary.pdf']),'-dpdf')
figure(AFig);
% print(fullfile(rc.fitOutputSummary, [date '_i613_i614_AuditoryCatchSummary.pdf']),'-dpdf')

