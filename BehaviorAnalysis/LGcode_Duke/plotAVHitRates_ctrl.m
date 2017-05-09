clear all
close all
earlyLT50 = 1;
lapseLT10 = 1;
% amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename_audCtrl, [], rc.indexTextCols);
% % av = behavParamsAV;
mice = unique(xd.Subject);
nMice = length(mice);

% fn = fullfile(rc.fitOutputSummary, [date '_CatchSummary.mat']);
fn = fullfile(rc.fitOutputSummary, ['11-Jan-2017_CatchSummary.mat']);
load(fn)

mnm = behavParamsMNM;

for imouse = 1:nMice;
    mouse_name = mnm(find(cat(1,mnm(:).subnum) == mice(imouse))).mouse;
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
    disp(size(use_ind))
    input = mouse(imouse).input(use_ind);
    input = concatenateStructures(input');
    
    
    missedIx = strcmp(input.trialOutcomeCell, 'ignore');
    successIx = strcmp(input.trialOutcomeCell, 'success');
    earliesIx = strcmp(input.trialOutcomeCell, 'failure');

    targetOri = chop(double(celleqel2mat_padded(input.tGratingDirectionDeg)),2);
    targetOris = unique(targetOri);
    
    cycTimeMs = unique(double(input.stimOffTimeMs)+double(input.stimOnTimeMs));
    
    
    nCyc = double(cell2mat(input.tCyclesOn));
    trialTimeMs = double(cycTimeMs*nCyc);
    lastStimTimeMs = double(cycTimeMs*(nCyc-1));
        
    targetReact = celleqel2mat_padded(input.reactTimeMs);    
    earlyReact = (celleqel2mat_padded(input.leverUpTimeMs) - celleqel2mat_padded(input.leverDownTimeMs)) - lastStimTimeMs;
    earlyReact(earlyReact < 100) = NaN;
    isBlock2 = cell2mat(input.isBlock2);
    
    % make cell array of data used for ea mouse
    tarRct{imouse} = targetReact;
    eRct{imouse} = earlyReact;
    
    tarDeg{imouse} = targetOri;
    
    sIx_all{imouse} = successIx;
    mIx_all{imouse} = missedIx;
    eIx_all{imouse} = earliesIx;
    bl2{imouse} = isBlock2;
    
    %bin trials by orientation
    [h_ori, bin_ori] = histc(targetOri, ori_edges);
    
    ori_bins = unique(bin_ori);
    avg_ori = zeros(1,max(ori_bins,[],2));
    sem_ori = zeros(1,max(ori_bins,[],2));
    n_ori = zeros(2,max(ori_bins,[],2));
    
    Hits_v = zeros(1,max(ori_bins,[],2));
    Misses_v = zeros(1,max(ori_bins,[],2));
    Earlies_v = zeros(1,max(ori_bins,[],2));
    rct_v = zeros(1,max(ori_bins,[],2));
    rctsem_v = zeros(1,max(ori_bins,[],2));
    Hits_av = zeros(1,max(ori_bins,[],2));
    Misses_av = zeros(1,max(ori_bins,[],2));
    Earlies_av = zeros(1,max(ori_bins,[],2));
    rct_av = zeros(1,max(ori_bins,[],2));
    rctsem_av = zeros(1,max(ori_bins,[],2));
    
    for ibin = 1:max(ori_bins,[],2)
        ind_v = intersect(find(bin_ori == ibin),find(isBlock2 == 0));
        avg_ori(ibin) = mean(targetOri(:,find(bin_ori == ibin)),2);
        sem_ori(ibin) = std(targetOri(:,find(bin_ori == ibin)),[],2)./sqrt(length(find(bin_ori == ibin)));
        if length(ind_v)>5
            Hits_v(ibin) = sum(successIx(ind_v),2);
            Misses_v(ibin) = sum(missedIx(ind_v),2);
            Earlies_v(ibin) = sum(earliesIx(ind_v),2);
            rct_v(ibin) = mean(targetReact(intersect(find(successIx),ind_v)),2);
            rctsem_v(ibin) = std(targetReact(intersect(find(successIx),ind_v)),[],2)./length(intersect(find(successIx),ind_v));
        end
        ind_av = intersect(find(bin_ori == ibin),find(isBlock2 == 1));
        if length(ind_av)>5
            Hits_av(ibin) = sum(successIx(ind_av),2);
            Misses_av(ibin) = sum(missedIx(ind_av),2);
            Earlies_av(ibin) = sum(earliesIx(ind_av),2);
            rct_av(ibin) = mean(targetReact(intersect(find(successIx),ind_av)),2);
            rctsem_av(ibin) = std(targetReact(intersect(find(successIx),ind_av)),[],2)./length(intersect(find(successIx),ind_av));
        end
    end
    [HR_v, ci95_HR_v] = binofit(Hits_v, Misses_v + Hits_v);
    [HR_av, ci95_HR_av] = binofit(Hits_av, Misses_av + Hits_av);
    [ER_v, ci95_ER_v] = binofit(Earlies_v, Earlies_v + Misses_v + Hits_v);
    [ER_av, ci95_ER_av] = binofit(Earlies_av, Earlies_av + Misses_av + Hits_av);
    
    n_ori = [Hits_v+Misses_v; Hits_av+Misses_av];
    
    bxSum = figure; 
    suptitle(['i' num2str(mouse_name) '- Vis and vis+tone trials']);
    
    subplot(2,2,1)
    h1 = errorbar(avg_ori, HR_v, HR_v - ci95_HR_v(:,1)', ci95_HR_v(:,2)' - HR_v);
    h1.Color = [0.5 0.5 0.5];
    h1.LineWidth = 3;
    h1.Marker = 'o';
    h1.MarkerFaceColor = [1 1 1];
    hold on 
    h2 = errorbar(avg_ori, HR_av, HR_av - ci95_HR_av(:,1)', ci95_HR_av(:,2)' - HR_av);
    h2.Color = [0 0 0];
    h2.LineWidth = 3;
    h2.Marker = 'o';
    h2.MarkerFaceColor = [1 1 1];
    set(gca, 'xscale', 'log');
    title('hit rate')
    xlabel('Orientation change (deg)')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([0 100])
    axis square
    legend(['V: n = ' num2str(sum(n_ori(1,:),2))], ['V+A: n = ' num2str(sum(n_ori(2,:),2))],'Location','Southeast');
%     print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_VisualCatchSummary.pdf']),'-dpdf')
 
    subplot(2,2,2) 
    h1 = errorbar(avg_ori, rct_v, rctsem_v);
    h1.Color = [0.5 0.5 0.5];
    h1.LineWidth = 3;
    h1.Marker = 'o';
    h1.MarkerFaceColor = [1 1 1];
    hold on 
    h2 = errorbar(avg_ori, rct_av, rctsem_av);
    h2.Color = [0 0 0];
    h2.LineWidth = 3;
    h2.Marker = 'o';
    h2.MarkerFaceColor = [1 1 1];
    set(gca, 'xscale', 'log');
    title('react time')
    xlabel('Orientation change (deg)')
    ylabel('RT (ms)')
%     ylim([250 350])
    xlim([0 100])
    axis square
    legend(['V: n = ' num2str(sum(n_ori(1,:),2))], ['V+A: n = ' num2str(sum(n_ori(2,:),2))],'Location','northeast');
%     print(fullfile(rc.fitOutputSummary, [date '_i' num2str(mouse_name) '_VisualCatchSummary.pdf']),'-dpdf')

%% plot HR by trial time
nTimes = 8;
trialTime_edges = [0:max(trialTimeMs)/nTimes:max(trialTimeMs)]+100;

[h_time, bin_time] = histc(trialTimeMs,trialTime_edges);

time_bins = unique(bin_time);
avg_trTime = zeros(1,length(time_bins));
hitsV_trTime = zeros(1,length(time_bins));
missesV_trTime = zeros(1,length(time_bins));
hitsAV_trTime = zeros(1,length(time_bins));
missesAV_trTime = zeros(1,length(time_bins));
rctV_trTime = zeros(1,length(time_bins));
rctV_sem_trTime = zeros(1,length(time_bins));
rctAV_trTime = zeros(1,length(time_bins));
rctAV_sem_trTime = zeros(1,length(time_bins));

for ibin = 1:length(time_bins)
    ind_v = intersect(find(bin_time == time_bins(ibin)), find(isBlock2 == 0));
    ind_av = intersect(find(bin_time == time_bins(ibin)), find(isBlock2 == 1));
    
    avg_trTime(ibin) = mean(trialTimeMs(find(bin_time == time_bins(ibin))));
    
    rctV_trTime(ibin) = mean(targetReact(intersect(find(successIx),ind_v)));
    rctV_sem_trTime(ibin) = nanstd(targetReact(intersect(find(successIx),ind_v)))/sqrt(length(intersect(find(successIx),ind_v)));
    rctAV_trTime(ibin) = mean(targetReact(intersect(find(successIx),ind_av)));
    rctAV_sem_trTime(ibin) = nanstd(targetReact(intersect(find(successIx),ind_av)))/sqrt(length(intersect(find(successIx),ind_av)));
    
    hitsV_trTime(ibin) = sum(successIx(ind_v));
    missesV_trTime(ibin) = sum(missedIx(ind_v));
    hitsAV_trTime(ibin) = sum(successIx(ind_av));
    missesAV_trTime(ibin) = sum(missedIx(ind_av));
end

    [HR_v_trTime, ci95_HR_v_trTime] = binofit(hitsV_trTime, missesV_trTime + hitsV_trTime);
    [HR_av_trTime, ci95_HR_av_trTime] = binofit(hitsAV_trTime, missesAV_trTime + hitsAV_trTime);
  
    %hit rate
    figure(bxSum)
    subplot(2,2,3)
    h1 = errorbar(floor(avg_trTime), HR_v_trTime, HR_v_trTime - ci95_HR_v_trTime(:,1)', ci95_HR_v_trTime(:,2)' - HR_v_trTime);
    h1.Color = [0.5 0.5 0.5];
    h1.LineWidth = 3;
    h1.LineStyle = 'none';
    h1.Marker = 'o';
    h1.MarkerFaceColor = [1 1 1];
    hold on 
    h2 = errorbar(floor(avg_trTime), HR_av_trTime, HR_av_trTime - ci95_HR_av_trTime(:,1)', ci95_HR_av_trTime(:,2)' - HR_av_trTime);
    h2.Color = [0 0 0];
    h2.LineWidth = 3;
    h2.LineStyle = 'none';
    h2.Marker = 'o';
    h2.MarkerFaceColor = [1 1 1];
    title('lasting effect of tone on hit rate')
    xlabel('time (ms)')
    ylabel('Hit rate')
    ylim([0 1])
    xlim([0 4000])
    axis square
    legend(['V: n = ' num2str(sum(n_ori(1,:),2))], ['V+A: n = ' num2str(sum(n_ori(2,:),2))],'Location','Southeast');
    
    %react time
    figure(bxSum)
    subplot(2,2,4)
    h1 = errorbar(floor(avg_trTime), rctV_trTime, rctV_sem_trTime);
    h1.Color = [0.5 0.5 0.5];
    h1.LineWidth = 3;
    h1.LineStyle = 'none';
    h1.Marker = 'o';
    h1.MarkerFaceColor = [1 1 1];
    hold on 
    h2  = errorbar(floor(avg_trTime), rctAV_trTime, rctAV_sem_trTime);
    h2.Color = [0 0 0];
    h2.LineWidth = 3;
    h2.LineStyle = 'none';
    h2.Marker = 'o';
    h2.MarkerFaceColor = [1 1 1];
    title('react time')
    xlabel('time (ms)')
    ylabel('RT (ms)')
%     ylim([250 350])
    xlim([0 4000])
    axis square
%     legend(['V: n = ' num2str(sum(n_ori(1,:),2))], ['V+A: n = ' num2str(sum(n_ori(2,:),2))],'Location','northeast');
    
    print(fullfile(rc.ashleyAnalysis, mouse_name,'behavior', 'bxSumAcrossDays'),'-dpdf','-fillpage')    

end