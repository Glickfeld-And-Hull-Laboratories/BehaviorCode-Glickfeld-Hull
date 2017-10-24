clear all
close all
rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
datasetStr = ['awFSAVdatasets_V1'];
fnout = fullfile(rc.caOutputDir,datasetStr,'behavior');
if ~exist(fnout,'dir')
    mkdir(fnout)
end

%%
minTrN_expt = 1;
minTrN_all = 10;
thresh = 30;
early_cutoff = 0.5;
lapse_cutoff = 0.9;
thresh_ms = 2000;
nTimes = 7;

nexp = xd.nRows;
mice = unique(xd.Subject);
nmice = length(mice);
colors_mice = parula(nmice);
% expCatch = logical(cell2mat({expt.catch}));


doExptStruct = 0;
if doExptStruct
    bxExptStruct_FSAV_train
else
    load(fullfile(fnout,'bxExpMat_training.mat'))
end

amp_edges = [0.003 0.007 0.02 0.05 0.14 0.37 1];
% amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

%%
ms = struct;
for im = 1:nmice
    ms(im).sn = mice(im);
end
visTargets_all = [];
invVisTargets_all = [];
invAudTargets_all = [];
audTargets_all = [];
successIx_all = [];
missedIx_all = [];
invHitIx_all = [];
invMissIx_all = [];
trLengthMs_all = [];
prevTrial_all = [];
prev2Trial_all = [];
prevVisTarget_all = [];

pctInv_expt = zeros(1,nexp);
invHR_expt = zeros(1,nexp);
valHR_expt = zeros(1,nexp);

doPlot = 0;
ntr1 = 0;
ntr2 = 0;
%% analyze behavior data from each experiment
hr_calc = 1;
for iexp = 1:nexp
    if isnan(bxExp(iexp).invType)
        continue
    else
        ntr1 = ntr1 + length(bxExp(iexp).sIx);
        if early_mat(iexp) > early_cutoff | HR_ori_mat(iexp) < lapse_cutoff | HR_amp_mat(iexp) < lapse_cutoff
            continue
        end
        ntr2 = ntr2 + length(bxExp(iexp).sIx);
        sn = bxExp(iexp).sn;
        sn_ind = find(mice == sn);
        cycLengthMs = double(bxExp(iexp).tOn + bxExp(iexp).tOff);

        visTargets = chop(double(bxExp(iexp).tVisTargets),2);
        audTargets = chop(double(bxExp(iexp).tAudTargets),2);

        nt = length(visTargets);
        
        trType = visTargets > 0;
        trType_shift = [NaN trType];
        prevTrial = trType_shift(1:length(trType));
        trType_shift = [NaN NaN trType];
        prev2Trial = trType_shift(1:length(trType));

        visTargets_shift = [NaN visTargets];
        prevVisTarget = visTargets_shift(1:length(visTargets));

        successIx = double(bxExp(iexp).sIx);
        missedIx = double(bxExp(iexp).mIx);
        invHitIx = double(bxExp(iexp).invHitIx);
        invMissIx = double(bxExp(iexp).invMissIx);
        nCyc = double(bxExp(iexp).trLength);
        trLengthMs = nCyc.*cycLengthMs;

        % sort vis catch and aud catch days, match invalid to valid target
        % types
        try
        if strcmp(bxExp(iexp).invType,'vis')
            invVisTargets = chop(double(bxExp(iexp).tInvTargets),2);
            invAudTargets = nan(1,nt);
            invTrials = invVisTargets > 0;
            invTars = round(unique(invVisTargets(invTrials)));
            matchValTrials = matchTrialsBx(invVisTargets,invTars,invTrials,visTargets);
            pctInv_expt(iexp) = sum(invTrials)./nt;
            invHR_expt(iexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
            valHR_expt(iexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));
            nInvTrials(iexp) = sum(invTrials);
        elseif strcmp(bxExp(iexp).invType,'aud')
            invVisTargets = nan(1,nt);
            invAudTargets = chop(double(bxExp(iexp).tInvTargets),2);  
            invTrials = invAudTargets > 0;
            invTars = unique(invAudTargets(invTrials));  
            matchValTrials = matchTrialsBx(invAudTargets,invTars,invTrials,audTargets); 
            pctInv_expt(iexp) = sum(invTrials)./nt;
            invHR_expt(iexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
            valHR_expt(iexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));  
            nInvTrials(iexp) = sum(invTrials);
        else
            invVisTargets = nan(1,nt);
            invAudTargets = nan(1,nt);
            pctInv_expt(iexp) = nan;
            invHR_expt(iexp) = nan;
            valHR_expt(iexp) = nan;
            nInvTrials(iexp) = nan;
        end
        catch
            invVisTargets = nan(1,nt);
            invAudTargets = nan(1,nt);
            pctInv_expt(iexp) = nan;
            invHR_expt(iexp) = nan;
            valHR_expt(iexp) = nan;
            nInvTrials(iexp) = nan;            
        end

        visTargets_all = cat(2,visTargets_all,visTargets);
        invVisTargets_all = cat(2,invVisTargets_all,invVisTargets);
        invAudTargets_all = cat(2,invAudTargets_all,invAudTargets);
        audTargets_all = cat(2,audTargets_all,audTargets);
        successIx_all = cat(2,successIx_all,successIx);
        missedIx_all = cat(2,missedIx_all,missedIx);
        invHitIx_all = cat(2,invHitIx_all,invHitIx);
        invMissIx_all = cat(2,invMissIx_all,invMissIx);
        trLengthMs_all = cat(2,trLengthMs_all,trLengthMs);
        prevTrial_all = cat(2,prevTrial_all,prevTrial);
        prev2Trial_all = cat(2,prev2Trial_all,prev2Trial);
        prevVisTarget_all = cat(2,prevVisTarget_all,prevVisTarget);

        if isfield(ms,'visTargets')           
            ms(sn_ind).visTargets = cat(2,ms(sn_ind).visTargets,visTargets);
            ms(sn_ind).invVisTargets = cat(2,ms(sn_ind).invVisTargets,invVisTargets);
            ms(sn_ind).invAudTargets = cat(2,ms(sn_ind).invAudTargets,invAudTargets);
            ms(sn_ind).audTargets = cat(2,ms(sn_ind).audTargets,audTargets);
            ms(sn_ind).successIx = cat(2,ms(sn_ind).successIx,successIx);
            ms(sn_ind).missedIx = cat(2,ms(sn_ind).missedIx,missedIx);
            ms(sn_ind).invHitIx = cat(2,ms(sn_ind).invHitIx,invHitIx);
            ms(sn_ind).invMissIx = cat(2,ms(sn_ind).invMissIx,invMissIx);
            ms(sn_ind).trLengthMs = cat(2,ms(sn_ind).trLengthMs,trLengthMs);
        else
            ms(sn_ind).visTargets = visTargets;
            ms(sn_ind).invVisTargets = invVisTargets;
            ms(sn_ind).invAudTargets = invAudTargets;
            ms(sn_ind).audTargets = audTargets;
            ms(sn_ind).successIx = successIx;
            ms(sn_ind).missedIx = missedIx;
            ms(sn_ind).invHitIx = invHitIx;
            ms(sn_ind).invMissIx = invMissIx;  
            ms(sn_ind).trLengthMs = trLengthMs;      
        end

        bxSumExptHR_FSAV_img
    end
end

%% plotting params
y_lim_hr = [0 1.1];
y_axis_hr = 0:0.2:1;
x_axis_ori = fliplr(chop(spoc(90,1,4),2));
x_tick_hr = [0 x_axis_ori];
x_tick_label_hr = [0 x_axis_ori];
x_lim_hr = [0 100];
%% behavior summary for each mouse

for im = 1:nmice
    
    % performance on visual trials across orientation
    [hr_v, ci95_v, avg_ori, sem_ori] = binHitRate(ms(im).visTargets,ms(im).successIx, ms(im).missedIx, ori_edges, minTrN_expt);
    [hr_av, ci95_av, avg_oric, sem_oric] = binHitRate(ms(im).invVisTargets,ms(im).invHitIx, ms(im).invMissIx, ori_edges, minTrN_expt);

    err_v = abs(bsxfun(@minus, ci95_v', hr_v));
    err_av = abs(bsxfun(@minus, ci95_av', hr_av));

    ind = ~isnan(hr_v);
    indc = ~isnan(hr_av);

    setFigParams4Print('landscape'); all_fig = figure;
    suptitle(sprintf('mouse: %s', num2str(ms(im).sn)));
    subplot 221
    h = errorbarxy(avg_ori(ind), hr_v(ind), sem_ori(ind), sem_ori(ind), err_v(1,ind), err_v(2,ind),{'ok', 'k', 'k'});
    h.hMain.MarkerFaceColor = 'k';
    h.hMain.LineStyle = '-';
    h.hMain.LineWidth = 1;
    hold on
    h = errorbarxy(avg_oric(indc), hr_av(indc), sem_oric(indc), sem_oric(indc), err_av(1,indc), err_av(2,indc),{'oc', 'c', 'c'});
    h.hMain.MarkerFaceColor = 'c';
    h.hMain.LineStyle = '-';
    h.hMain.LineWidth = 1;
    h.hMain.Parent.XScale = 'log';
    figXAxis([],'Orientation change (deg)',x_lim_hr,x_tick_hr,x_tick_label_hr);
    figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
    figAxForm([])
    title('visual trials')
    
    % performance on auditory trials across orientation
    [hr_a, ci95_a, avg_amp, sem_amp] = binHitRate(ms(im).audTargets,ms(im).successIx, ms(im).missedIx, amp_edges, minTrN_expt);
    if sum(ms(im).invAudTargets > 0) == 0
        [hr_va, ci95_va, avg_ampc, sem_ampc] = deal(nan);
    else
        [hr_va, ci95_va, avg_ampc, sem_ampc] = binHitRate(ms(im).invAudTargets,ms(im).invHitIx, ms(im).invMissIx, amp_edges, minTrN_expt);
    end
    err_a = abs(bsxfun(@minus, ci95_a', hr_a));
    err_va = abs(bsxfun(@minus, ci95_va', hr_va));
    
    ind = ~isnan(hr_a);
    indc = ~isnan(hr_va);

    subplot 222
    h = errorbarxy(avg_amp(ind), hr_a(ind), sem_amp(ind), sem_amp(ind), err_a(1,ind), err_a(2,ind),{'ok', 'k', 'k'});
    h.hMain.MarkerFaceColor = 'k';
    h.hMain.LineStyle = '-';
    h.hMain.LineWidth = 1;
    hold on
    if sum(ms(im).invAudTargets > 0) ~= 0
        h = errorbarxy(avg_ampc(indc), hr_va(indc), sem_ampc(indc), sem_ampc(indc), err_va(1,indc), err_va(2,indc),{'oc', 'c', 'c'});
        h.hMain.MarkerFaceColor = 'c';
        h.hMain.LineStyle = '-';
        h.hMain.LineWidth = 1;
    end
    h.hMain.Parent.XScale = 'log';
    figXAxis([],'volume (% of max)',[0.001 0.5], [0.001,0.01,0.1,0.5],[0.001,0.01,0.1,0.5]);
    figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
    figAxForm([])
    title('auditory trials')
        
    % performance on visual trials across trial length
    trL = ms(im).trLengthMs;
    trialTime_edges = [0:max(trL)/nTimes:max(trL)]+100;

    [hr_v, ci95_v, avg_time, sem_time] = binHitRate(trL,ms(im).successIx, ms(im).missedIx, trialTime_edges, minTrN_expt);
    [v_fit,v_trL] = lgtFitHR(trL',ms(im).visTargets' > 0,ms(im).successIx',ms(im).missedIx');
    [v_trL,sort_ind] = sort(v_trL);
    v_fit = v_fit(sort_ind);
    [hr_av, ci95_av, avg_timec, sem_timec] = binHitRate(trL,ms(im).invHitIx, ms(im).invMissIx, trialTime_edges, minTrN_expt);
    [av_fit,av_trL] = lgtFitHR(trL',ms(im).invVisTargets' > 0,ms(im).invHitIx',ms(im).invMissIx');
    [av_trL,sort_ind] = sort(av_trL);
    av_fit = av_fit(sort_ind);
    err_v = abs(bsxfun(@minus, ci95_v', hr_v));
    err_av = abs(bsxfun(@minus, ci95_av', hr_av));

    ind = ~isnan(hr_v);
    indc = ~isnan(hr_av);
    x_tick = unique([0 chop(avg_time,2)]);
    x_tick_label = unique([0 chop(avg_time./1000,2)]);
    
    subplot 223
    h = errorbarxy(avg_time(ind), hr_v(ind), sem_time(ind), sem_time(ind), err_v(1,ind), err_v(2,ind),{'ok', 'k', 'k'});
    h.hMain.MarkerFaceColor = 'k';
%     h.hMain.LineStyle = '-';
%     h.hMain.LineWidth = 1;
    hold on
    h = plot(v_trL,v_fit,'k-');
    h.LineWidth = 1;
    hold on
    h = errorbarxy(avg_timec(indc), hr_av(indc), sem_timec(indc), sem_timec(indc), err_av(1,indc), err_av(2,indc),{'oc', 'c', 'c'});
    h.hMain.MarkerFaceColor = 'c';
    hold on
    h = plot(av_trL,av_fit,'c-');
    h.LineWidth = 1;
%     h.hMain.LineStyle = '-';
%     h.hMain.LineWidth = 1;
    figXAxis([],'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
    figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
    figAxForm([])
    title('visual trials')
    
    % performance on visual trials >1.5s 
    time_ind = trL > thresh_ms;
    [hr_v, ci95_v, avg_ori, sem_ori] = binHitRate(ms(im).visTargets(time_ind),ms(im).successIx(time_ind), ms(im).missedIx(time_ind), ori_edges, minTrN_expt);
    [hr_av, ci95_av, avg_oric, sem_oric] = binHitRate(ms(im).invVisTargets(time_ind),ms(im).invHitIx(time_ind), ms(im).invMissIx(time_ind), ori_edges, minTrN_expt);

    err_v = abs(bsxfun(@minus, ci95_v', hr_v));
    err_av = abs(bsxfun(@minus, ci95_av', hr_av));

    ind = ~isnan(hr_v);
    indc = ~isnan(hr_av);
    
    subplot 224
    h = errorbarxy(avg_ori(ind), hr_v(ind), sem_ori(ind), sem_ori(ind), err_v(1,ind), err_v(2,ind),{'ok', 'k', 'k'});
    h.hMain.MarkerFaceColor = 'k';
    h.hMain.LineStyle = '-';
    h.hMain.LineWidth = 1;
    hold on
    if sum(indc) ~= 0
        h = errorbarxy(avg_oric(indc), hr_av(indc), sem_oric(indc), sem_oric(indc), err_av(1,indc), err_av(2,indc),{'oc', 'c', 'c'});
        h.hMain.MarkerFaceColor = 'c';
        h.hMain.LineStyle = '-';
        h.hMain.LineWidth = 1;
    end
    h.hMain.Parent.XScale = 'log';
    figXAxis([],'Orientation change (deg)',x_lim_hr,x_tick_hr,x_tick_label_hr);
    figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
    figAxForm([])
    title(sprintf('trial length > %s s',num2str(thresh_ms/1000)))
        
    print(fullfile(fnout,sprintf('ms%s',num2str(ms(im).sn))),'-dpdf','-fillpage')
end

%% analyze behavior data across all experiments
hr_calc = 2;

bxSumExptHR_FSAV_img

ori_ind = ~isnan(HR_v);
oric_ind = ~isnan(HR_av);

% hit rate - valid and invalid visual trials
setFigParams4Print('portrait'); all_fig = figure;
suptitle('all imaging expt')
figure(all_fig);
sp2 = subplot(2,2,1);
h = errorbarxy(avg_ori(ori_ind), HR_v(ori_ind), sem_ori(ori_ind), sem_ori(ori_ind), HR_v(ori_ind) - ci95_HR_v((ori_ind),1)', ci95_HR_v((ori_ind),2)' - HR_v(ori_ind), {'ok', 'k', 'k'});
h.hMain.MarkerFaceColor  = 'k';
h.hMain.LineStyle = '-';
h.hMain.LineWidth = 3;
hold on
h = errorbarxy(avg_oric(oric_ind), HR_av(oric_ind), sem_oric(oric_ind), sem_oric(oric_ind), HR_av(oric_ind) - ci95_HR_av((oric_ind),1)', ci95_HR_av((oric_ind),2)' - HR_av(oric_ind), {'oc', 'c', 'c'});
h.hMain.MarkerFaceColor  = 'c';
h.hMain.LineStyle = '-';
h.hMain.LineWidth = 3;
% set(gca, 'xscale', 'log');
y_tick = y_axis_hr;
figXAxis(h.hMain.Parent,'Orientation change (deg)',x_lim_hr,x_tick_hr,x_tick_label_hr);
figYAxis(h.hMain.Parent,'HR',[0 1.1],y_tick);
figAxForm(sp2)
% h.hMain.Parent.XTick = [ori_edges(3) ori_edges(end)];
% h.hMain.Parent.XTickLabel = ;
% xlabel('Orientation change (deg)')
% ylabel('Hit rate')
% ylim([0 1.1])
% xlim([10 100])
title({['val = ' num2str(sum(n_ori(1,:),2))], ['inv = ' num2str(sum(n_ori(2,:),2))]});
% axis square

% auditory trials
[hr_a, ci95_a, avg_amp, sem_amp] = binHitRate(audTargets_all,successIx_all,missedIx_all, amp_edges, minTrN_all);
[hr_va, ci95_va, avg_ampc, sem_ampc] = binHitRate(invAudTargets_all, invHitIx_all, invMissIx_all,amp_edges, minTrN_all);

err_a = abs(bsxfun(@minus, ci95_a', hr_a));
err_va = abs(bsxfun(@minus, ci95_va', hr_va));

ind = ~isnan(hr_a);
indc = ~isnan(hr_va);

figure(all_fig);
sp2 = subplot(2,2,2);
h = errorbarxy(avg_amp(ind), hr_a(ind), sem_amp(ind), sem_amp(ind), err_a(1,ind), err_a(1,ind), {'ok', 'k', 'k'});
h.hMain.MarkerFaceColor  = 'k';
h.hMain.LineStyle = '-';
h.hMain.LineWidth = 3;
hold on
h = errorbarxy(avg_ampc(indc), hr_va(indc), sem_ampc(indc), sem_ampc(indc), err_va(1,indc), err_va(1,indc), {'oc', 'c', 'c'});
h.hMain.MarkerFaceColor  = 'c';
h.hMain.LineStyle = '-';
h.hMain.LineWidth = 3;
h.hMain.Parent.XScale = 'log';
figXAxis([],'volume (% of max)',[0.001 1], [0.001,0.01,0.1,1],[0.001,0.01,0.1,1]);
figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
figAxForm([])
title('auditory trials')


% hit rate by trial length
bxSumTiming_FSAV_img

[h_vistime, p_vistime, chistat_vistime] = binChiHR(hits_trTime,misses_trTime,invHit_trTime,invMiss_trTime);
h_vistime = logical(h_vistime);

figure(all_fig)
sp3 = subplot(2,2,3);
h = errorbar(avg_trTime,HR_trTime,HR_trTime-ci_95_HR_trTime(:,1)',ci_95_HR_trTime(:,2)'-HR_trTime,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(avg_trTime,FR_trTime,FR_trTime-ci_95_FR_trTime(:,1)',ci_95_FR_trTime(:,2)'-FR_trTime,'co');
h.MarkerFaceColor = 'c';
hold on
h = plot(avg_trTime(h_vistime),FR_trTime(h_vistime),'co');
h.MarkerFaceColor = 'c';
h.MarkerEdgeColor = 'k';
hold on

x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[trialTime_edges(1) trialTime_edges(end)],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(sp3)
title('visual trials')

% timing & auditory trials
trL = trLengthMs_all;
trialTime_edges = [0:max(trL)/nTimes:max(trL)]+100;

[hr_a, ci95_a, avg_time, sem_time,hits_a,misses_a] = binHitRate(trL,successIx_all,missedIx_all, trialTime_edges, minTrN_all);
% [v_fit,v_trL] = lgtFitHR(trL',ms(im).visTargets' > 0,ms(im).successIx',ms(im).missedIx');
% [v_trL,sort_ind] = sort(v_trL);
% v_fit = v_fit(sort_ind);
[hr_va, ci95_va, avg_timec, sem_timec,hits_va,misses_va] = binHitRate(trL,invHitIx_all, invMissIx_all, trialTime_edges, minTrN_all);
% [av_fit,av_trL] = lgtFitHR(trL',ms(im).invVisTargets' > 0,ms(im).invHitIx',ms(im).invMissIx');
% [av_trL,sort_ind] = sort(av_trL);
% av_fit = av_fit(sort_ind);
[h_audtime, p_audtime, chistat_audtime] = binChiHR(hits_a,misses_a,hits_va,misses_va);
h_audtime = logical(h_audtime);
err_v = abs(bsxfun(@minus, ci95_a', hr_a));
err_av = abs(bsxfun(@minus, ci95_va', hr_va));

ind = ~isnan(hr_a);
indc = ~isnan(hr_va);

subplot 224
h = errorbarxy(avg_time(ind), hr_a(ind), sem_time(ind), sem_time(ind), err_v(1,ind), err_v(2,ind),{'ok', 'k', 'k'});
h.hMain.MarkerFaceColor = 'k';
% h.hMain.LineStyle = '-';
% h.hMain.LineWidth = 1;
hold on
% h = plot(v_trL,v_fit,'k-');
% h.LineWidth = 1;
% hold on
h = errorbarxy(avg_timec(indc), hr_va(indc), sem_timec(indc), sem_timec(indc), err_av(1,indc), err_av(2,indc),{'oc', 'c', 'c'});
h.hMain.MarkerFaceColor = 'c';
hold on
h = plot(avg_timec(h_audtime),hr_va(h_audtime),'co');
h.MarkerFaceColor = 'c';
h.MarkerEdgeColor = 'k';
% h = plot(av_trL,av_fit,'c-');
% h.LineWidth = 1;
% h.hMain.LineStyle = '-';
% h.hMain.LineWidth = 1;
figXAxis([],'trial length (s)',[trialTime_edges(1) trialTime_edges(end)],x_tick,x_tick_label);
figYAxis([],'Hit Rate (%)',y_lim_hr,y_axis_hr,y_axis_hr*100);
figAxForm([])
title('auditory trials')

figure(all_fig)
print(fullfile(fnout,'all_expt_train'),'-dpdf','-fillpage')

%% behavior dep on previous trial
% hit rate - vis trials sorted by preceeding trial type
ori_ind = ~isnan(HR_oripv);
oric_ind = ~isnan(FR_oricpv);
setFigParams4Print('landscape'); figure
suptitle('visual trial performance, dependence on trial history')
subplot 121
h = errorbar(avg_ori(ori_ind), HR_oripv(ori_ind), HR_oripv(ori_ind)-ci95_HR_oripv(ori_ind,1)', HR_oripv(ori_ind)-ci95_HR_oripv(ori_ind,2)','ko');
h.MarkerFaceColor  = 'k';
h.LineStyle = '-';
h.LineWidth = 1;
leg(1) = h;
hold on
h = errorbar(avg_ori(ori_ind), HR_oripa(ori_ind), HR_oripa(ori_ind)-ci95_HR_oripa(ori_ind,1)', HR_oripa(ori_ind)-ci95_HR_oripa(ori_ind,2)','ko');
h.MarkerFaceColor  = [0.5 0.5 0.5];
h.Color = [0.5 0.5 0.5];
h.LineStyle = '-';
h.LineWidth = 1;
leg(2) = h;
hold on
h = errorbar(avg_oric(oric_ind), FR_oricpv(oric_ind), FR_oricpv(oric_ind)-ci95_FR_oricpv(oric_ind,1)', FR_oricpv(oric_ind)-ci95_FR_oricpv(oric_ind,2)','co');
h.MarkerFaceColor  = 0.5*h.Color;
h.Color = 0.5*h.Color;
h.LineStyle = '-';
h.LineWidth = 1;
leg(3) = h;
h = errorbar(avg_oric(oric_ind), FR_oricpa(oric_ind), FR_oricpa(oric_ind)-ci95_FR_oricpa(oric_ind,1)', FR_oricpa(oric_ind)-ci95_FR_oricpa(oric_ind,2)','co');
h.MarkerFaceColor  = h.Color;
h.LineStyle = '-';
h.LineWidth = 1;
leg(4) = h;
h.Parent.XScale = 'log';
x_tick = [0 x_axis_ori];
x_tick_label = [0 x_axis_ori];
y_tick = y_axis_hr;
figXAxis(h.Parent,'Orientation change (deg)',[0 100],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_tick);
figAxForm(h.Parent)
legend(leg,{'vis->vis';'aud->vis';'vis->inv vis';'aud->inv vis'},'location','northwest')

subplot 122
h = errorbar(avg_trTime,HR_pv_trTime,HR_pv_trTime-ci_95_HR_pv_trTime(:,1)',ci_95_HR_pv_trTime(:,2)'-HR_pv_trTime,'ko');
% h = plot(avg_trTime,HR_pv_trTime,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(avg_trTime,HR_pa_trTime,HR_pa_trTime-ci_95_HR_pa_trTime(:,1)',ci_95_HR_pa_trTime(:,2)'-HR_pa_trTime,'ko');
% h = plot(avg_trTime,HR_pa_trTime,'ko');
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h = errorbar(avg_trTime,FR_pv_trTime,FR_pv_trTime-ci_95_FR_pv_trTime(:,1)',ci_95_FR_pv_trTime(:,2)'-FR_pv_trTime,'co');
% h = plot(avg_trTime,FR_pv_trTime,'co');
h.MarkerFaceColor = 0.5*h.Color;
h.Color = 0.5*h.Color;
h = errorbar(avg_trTime,FR_pa_trTime,FR_pa_trTime-ci_95_FR_pa_trTime(:,1)',ci_95_FR_pa_trTime(:,2)'-FR_pa_trTime,'co');
% h = plot(avg_trTime,FR_pa_trTime,'co');
h.MarkerFaceColor = 'c';
hold on
% h.Parent.XScale = 'log';
x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(h.Parent)
print(fullfile(fnout,'all_expt_train_depPrevTrial'),'-dpdf','-fillpage')

%% behavior dep on previous 2 trials
% hit rate - vis trials sorted by preceeding trial type
ori_ind = ~isnan(HR_oripvv);
oric_ind = ~isnan(FR_oricpvv);
setFigParams4Print('landscape'); figure
suptitle('visual trial performance, dependence on trial history')
subplot 121
h = errorbar(avg_ori(ori_ind), HR_oripvv(ori_ind), HR_oripvv(ori_ind)-ci95_HR_oripvv(ori_ind,1)', HR_oripvv(ori_ind)-ci95_HR_oripvv(ori_ind,2)','ko');
h.MarkerFaceColor  = 'k';
h.LineStyle = '-';
h.LineWidth = 1;
leg(1) = h;
hold on
h = errorbar(avg_ori(ori_ind), HR_oripaa(ori_ind), HR_oripaa(ori_ind)-ci95_HR_oripaa(ori_ind,1)', HR_oripaa(ori_ind)-ci95_HR_oripaa(ori_ind,2)','ko');
h.MarkerFaceColor  = [0.5 0.5 0.5];
h.Color = [0.5 0.5 0.5];
h.LineStyle = '-';
h.LineWidth = 1;
leg(2) = h;
hold on
h = errorbar(avg_oric(oric_ind), FR_oricpvv(oric_ind), FR_oricpvv(oric_ind)-ci95_FR_oricpvv(oric_ind,1)', FR_oricpvv(oric_ind)-ci95_FR_oricpvv(oric_ind,2)','co');
h.MarkerFaceColor  = 0.5*h.Color;
h.Color = 0.5*h.Color;
h.LineStyle = '-';
h.LineWidth = 1;
leg(3) = h;
h = errorbar(avg_oric(oric_ind), FR_oricpaa(oric_ind), FR_oricpaa(oric_ind)-ci95_FR_oricpaa(oric_ind,1)', FR_oricpaa(oric_ind)-ci95_FR_oricpaa(oric_ind,2)','co');
h.MarkerFaceColor  = h.Color;
h.LineStyle = '-';
h.LineWidth = 1;
leg(4) = h;
h.Parent.XScale = 'log';
x_tick = [0 x_axis_ori];
x_tick_label = [0 x_axis_ori];
y_tick = y_axis_hr;
figXAxis(h.Parent,'Orientation change (deg)',[0 100],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_tick);
figAxForm(h.Parent)
legend(leg,{'vis->vis->vis';'aud->aud->vis';'vis->vis->inv vis';'aud->aud->inv vis'},'location','northwest')

subplot 122
h = errorbar(avg_trTime,HR_pvv_trTime,HR_pvv_trTime-ci_95_HR_pvv_trTime(:,1)',ci_95_HR_pvv_trTime(:,2)'-HR_pvv_trTime,'ko');
% h = plot(avg_trTime,HR_pv_trTime,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(avg_trTime,HR_paa_trTime,HR_paa_trTime-ci_95_HR_paa_trTime(:,1)',ci_95_HR_paa_trTime(:,2)'-HR_paa_trTime,'ko');
% h = plot(avg_trTime,HR_pa_trTime,'ko');
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
h = errorbar(avg_trTime,FR_pvv_trTime,FR_pvv_trTime-ci_95_FR_pvv_trTime(:,1)',ci_95_FR_pvv_trTime(:,2)'-FR_pvv_trTime,'co');
% h = plot(avg_trTime,FR_pv_trTime,'co');
h.MarkerFaceColor = 0.5*h.Color;
h.Color = 0.5*h.Color;
h = errorbar(avg_trTime,FR_paa_trTime,FR_paa_trTime-ci_95_FR_paa_trTime(:,1)',ci_95_FR_paa_trTime(:,2)'-FR_paa_trTime,'co');
% h = plot(avg_trTime,FR_pa_trTime,'co');
h.MarkerFaceColor = 'c';
hold on
% h.Parent.XScale = 'log';
x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(h.Parent)
print(fullfile(fnout,'all_expt_train_depPrev2Trial'),'-dpdf','-fillpage')