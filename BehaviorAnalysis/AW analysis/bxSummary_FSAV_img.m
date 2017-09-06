clear all
close all
ds = '_V1som';
%%
minTrN_expt = 1;
minTrN_all = 1;
nTimes = 5;

datasetStr = ['awFSAVdatasets' ds];
eval(datasetStr)
rc = behavConstsAV;
fnout = fullfile(rc.caOutputDir,datasetStr,'behavior');
if ~exist(fnout,'dir')
    mkdir(fnout)
end
nexp = size(expt,2);
mice = unique({expt.SubNum});
nmice = length(mice);
colors_mice = parula(nmice);
expCatch = logical(cell2mat({expt.catch}));


doExptStruct = 1;
if doExptStruct
    bxExptStruct_FSAV_img
end

amp_edges = [0.003 0.007 0.02 0.05 0.14 0.37 1];
% amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

%%
ms_str = struct;
for im = 1:nmice
    ms_str(im).sn = mice(im);
end
visTargets_all = [];
invVisTargets_all = [];
audTargets_all = [];
successIx_all = [];
missedIx_all = [];
invHitIx_all = [];
invMissIx_all = [];
trLengthMs_all = [];
prevTrial_all = [];
prev2Trial_all = [];
prevVisTarget_all = [];

invHR_expt = zeros(1,nexp);
valHR_expt = zeros(1,nexp);

doPlot = 0;
%% analyze behavior data from each experiment
hr_calc = 1;
for iexp = 1:nexp
    if expt(iexp).catch
    title_str = [expt(iexp).SubNum '-' expt(iexp).date];
    sn = bxExp(iexp).sn;
    sn_ind = find(strcmp(mice,sn));
    cycLengthMs = double(bxExp(iexp).tOn + bxExp(iexp).tOff);
    
    visTargets = chop(double(bxExp(iexp).tVisTargets),2);
    invVisTargets = chop(double(bxExp(iexp).tInvVisTargets),2);
    audTargets = chop(double(bxExp(iexp).tAudTargets),2);
    
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
    
    invTrials = invVisTargets > 0;
    invTars = round(unique(invVisTargets(invTrials)));
    
    
    matchValTrials = matchTrialsBx(invVisTargets,invTars,invTrials,visTargets);
    
    invHR_expt(iexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
    valHR_expt(iexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));
    
    visTargets_all = cat(2,visTargets_all,visTargets);
    invVisTargets_all = cat(2,invVisTargets_all,invVisTargets);
    audTargets_all = cat(2,audTargets_all,audTargets);
    successIx_all = cat(2,successIx_all,successIx);
    missedIx_all = cat(2,missedIx_all,missedIx);
    invHitIx_all = cat(2,invHitIx_all,invHitIx);
    invMissIx_all = cat(2,invMissIx_all,invMissIx);
    trLengthMs_all = cat(2,trLengthMs_all,trLengthMs);
    prevTrial_all = cat(2,prevTrial_all,prevTrial);
    prev2Trial_all = cat(2,prev2Trial_all,prev2Trial);
    prevVisTarget_all = cat(2,prevVisTarget_all,prevVisTarget);
    
    if isfield(ms_str,'visTargets')
        if isempty(ms_str(sn_ind).visTargets)
            ms_str(sn_ind).visTargets = visTargets;
            ms_str(sn_ind).invVisTargets = invVisTargets;
            ms_str(sn_ind).audTargets = audTargets;
            ms_str(sn_ind).successIx = successIx;
            ms_str(sn_ind).missedIx = missedIx;
            ms_str(sn_ind).invHitIx = invHitIx;
            ms_str(sn_ind).invMissIx = invMissIx;
            ms_str(sn_ind).trLengthMs = trLengthMs;
        else            
            ms_str(sn_ind).visTargets = cat(2,ms_str(sn_ind).visTargets,visTargets);
            ms_str(sn_ind).invVisTargets = cat(2,ms_str(sn_ind).invVisTargets,invVisTargets);
            ms_str(sn_ind).audTargets = cat(2,ms_str(sn_ind).audTargets,audTargets);
            ms_str(sn_ind).successIx = cat(2,ms_str(sn_ind).successIx,successIx);
            ms_str(sn_ind).missedIx = cat(2,ms_str(sn_ind).missedIx,missedIx);
            ms_str(sn_ind).invHitIx = cat(2,ms_str(sn_ind).invHitIx,invHitIx);
            ms_str(sn_ind).invMissIx = cat(2,ms_str(sn_ind).invMissIx,invMissIx);
            ms_str(sn_ind).trLengthMs = cat(2,ms_str(sn_ind).trLengthMs,trLengthMs);
        end
    else
            ms_str(sn_ind).visTargets = visTargets;
            ms_str(sn_ind).invVisTargets = invVisTargets;
            ms_str(sn_ind).audTargets = audTargets;
            ms_str(sn_ind).successIx = successIx;
            ms_str(sn_ind).missedIx = missedIx;
            ms_str(sn_ind).invHitIx = invHitIx;
            ms_str(sn_ind).invMissIx = invMissIx;  
            ms_str(sn_ind).trLengthMs = trLengthMs;      
    end
    
    bxSumExptHR_FSAV_img
    bxSumExptHR_audtr_FSAV_img
    
%     % sound amplitude values are currently incorrect
    if doPlot
        expt_fig = figure; setFigParams4Print('portrait')
        suptitle(title_str)
        subplot(4,2,1)
        errorbarxy(avg_ori, HR_ori, sem_ori, sem_ori, HR_ori - ci95_HR_v(:,1)', ci95_HR_v(:,2)' - HR_ori, {'ok', 'k', 'k'});
        hold on
        errorbarxy(avg_oric, FR_ori, sem_oric, sem_oric, FR_ori - ci95_HR_av(:,1)', ci95_HR_av(:,2)' - FR_ori, {'oc', 'c', 'c'});
        set(gca, 'xscale', 'log');
        xlabel('Orientation change (deg)')
        ylabel('Hit rate')
        ylim([0 1])
        xlim([10 100])
        title({['val = ' num2str(sum(n_ori(1,:),2))], ['inv = ' num2str(sum(n_ori(2,:),2))]});
        axis square
        
        figure(expt_fig);
        subplot(4,2,2)
        errorbarxy(avg_amp, HR_amp, sem_amp, sem_amp, HR_amp - ci95_HR_amp(:,1)', ci95_HR_amp(:,2)' - HR_amp, {'ok', 'k', 'k'});
        hold on
%         errorbarxy(avg_oric, FR_ori, sem_oric, sem_oric, FR_ori - ci95_HR_av(:,1)', ci95_HR_av(:,2)' - FR_ori, {'oc', 'c', 'c'});
        set(gca, 'xscale', 'log');
        xlabel('Amplitude (% max)')
        ylabel('Hit rate')
        ylim([0 1])
%         xlim([0 0.2])
        title({['val = ' num2str(sum(n_amp(1,:),2))]});
        axis square
       
    if iexp == 14;
        exMsBxFig = figure;
        suptitle(title_str)
        sp1 = subplot(1,2,1);
        h = plot(avg_ori, HR_ori, 'ko');
        h.MarkerFaceColor = 'k';
        hold on
        h = plot(avg_oric, FR_ori, 'co');
        h.MarkerFaceColor = 'c';
        xlabel('Orientation change (deg)')
        ylabel('Hit rate')
        x_axis_ori = fliplr(chop(spoc(90,1,4),2));
        sp1.XTick = [0 x_axis_ori];
        sp1.XTickLabel = [0 x_axis_ori];
%         sp1.XScale = 'log';
        ylim([0 1])
        xlim([0 100])
        title({['val = ' num2str(sum(n_ori(1,:),2))], ['inv = ' num2str(sum(n_ori(2,:),2))]});
        axis square
        
        figure(exMsBxFig);
        sp2 = subplot(1,2,2);
        h = plot(avg_amp, HR_amp, 'ko');
        hold on
%         errorbarxy(avg_oric, FR_ori, sem_oric, sem_oric, FR_ori - ci95_HR_av(:,1)', ci95_HR_av(:,2)' - FR_ori, {'oc', 'c', 'c'});
%         set(gca, 'xscale', 'log');
        x_axis_amp = fliplr(chop(spoc(0.2,1,4),2));
        sp2.XTick = [0 x_axis_amp];
        sp2.XTickLabel = [0 x_axis_amp];
        xlabel('Amplitude (% max)')
        ylabel('Hit rate')
        ylim([0 1])
%         xlim([0 0.2])
        title({['val = ' num2str(sum(n_amp(1,:),2))]});
        axis square
        print([fnout '\exDay'],'-dpdf','-fillpage')
    end

        figure(expt_fig);
        subplot(4,2,7)
        x_run = 1:expt(iexp).nrun+1;
        far = cat(2,bxExp(iexp).far_run,early_mat(iexp));
        h = plot(x_run,far,'ko');
        h.MarkerFaceColor = 'k';
        xlim([0 x_run(end)+1])
        ylim([0 1])
        h.Parent.XTick = x_run;
        h.Parent.XTickLabel = cat(2,strsplit(num2str(x_run(1:end-1))), 'all');
        xlabel('expt run')
        ylabel('% early release')
        axis square

        figure(expt_fig);
        subplot(4,2,8)
        hrve = cat(2,bxExp(iexp).hr_run,HR_ori_mat(iexp));
        h = plot(x_run,hrve,'ko');
        h.MarkerFaceColor = 'k';
        hold on
        hrae = HR_amp_mat(iexp);
        plot(x_run(end),hrae,'ko');
        xlim([0 x_run(end)+1])
        ylim([0 1])
        h.Parent.XTick = x_run;
        h.Parent.XTickLabel = cat(2,strsplit(num2str(x_run(1:end-1))), 'all');
        xlabel('expt run')
        ylabel('hit rate easy trials')
        axis square
        figure(expt_fig)
        print(fullfile(fnout,title_str),'-dpdf','-fillpage')
    end
    end
end

y_axis_hr = 0:0.2:1;
all_fig = figure; setFigParams4Print('portrait')
suptitle('all imaging expt')
sp1 = subplot(2,2,1);
for iexp = 1:nexp
    
    sn_ind = find(strcmp(mice,bxExp(iexp).sn));
    if expt(iexp).catch
        l = [valHR_expt(iexp)' invHR_expt(iexp)'];
        h =  plot(1:2,l,'o-');
        h.Color = [0.5 0.5 0.5];
%         if ~expt(iexp).catchRew
%             h.MarkerEdgeColor = [1 0 0];
%         else
%             h.MarkerEdgeColor = [1 1 1];
%         end
        h.MarkerFaceColor = colors_mice(sn_ind,:);
%         h.MarkerFaceColor = colors_exp(iexp,:);
    end
    hold on
end

title('HR each exp')
legend({expt(expCatch).SubNum},'location','eastoutside')

figXAxis(sp1,[],[0.75 2.25],[1 2],{'valid','invalid'});
figYAxis(sp1,'match HR',[0 1.1],y_axis_hr);
sp1.TickDir = 'out';
sp1.Box = 'off';
%% analyze behavior data across all experiments
hr_calc = 2;

bxSumExptHR_FSAV_img

ori_ind = ~isnan(HR_v);
oric_ind = ~isnan(HR_av);
x_axis_ori = fliplr(chop(spoc(90,1,4),2));

% hit rate - valid and invalid visual trials
figure(all_fig);
sp2 = subplot(2,2,2);
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
x_tick = [0 x_axis_ori];
x_tick_label = [0 x_axis_ori];
y_tick = y_axis_hr;
figXAxis(h.hMain.Parent,'Orientation change (deg)',[0 100],x_tick,x_tick_label);
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

% hit rate by trial length
bxSumTiming_FSAV_img

figure(all_fig)
sp3 = subplot(2,2,3);
h = errorbar(avg_trTime,HR_trTime,HR_trTime-ci_95_HR_trTime(:,1)',ci_95_HR_trTime(:,2)'-HR_trTime,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(avg_trTime,FR_trTime,FR_trTime-ci_95_FR_trTime(:,1)',ci_95_FR_trTime(:,2)'-FR_trTime,'co');
h.MarkerFaceColor = 'c';
hold on

x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(sp3)
title('visual trials')

figure(all_fig)
print(fullfile(fnout,'all_expt_img'),'-dpdf','-fillpage')

%% behavior dep on previous trial
% hit rate - vis trials sorted by preceeding trial type

figure; setFigParams4Print('landscape')
subplot 121
% h = errorbarxy(avg_ori(ori_ind), HR_oripv(ori_ind), sem_ori(ori_ind), sem_ori(ori_ind), HR_oripv(ori_ind) - ci95_HR_oripv((ori_ind),1)', ci95_HR_oripv((ori_ind),2)' - HR_oripv(ori_ind), {'ok', 'k', 'k'});
h = plot(avg_ori(ori_ind), HR_oripvv(ori_ind),'ko');
h.MarkerFaceColor  = 'k';
h.LineStyle = '-';
h.LineWidth = 2;
leg(1) = h;
hold on
% h = errorbarxy(avg_ori(ori_ind), HR_oripa(ori_ind), sem_ori(ori_ind), sem_ori(ori_ind), HR_oripa(ori_ind) - ci95_HR_oripa((ori_ind),1)', ci95_HR_oripa((ori_ind),2)' - HR_oripa(ori_ind), {'ok', 'k', 'k'});
h = plot(avg_ori(ori_ind), HR_oripaa(ori_ind),'ko');
h.MarkerFaceColor  = [0.5 0.5 0.5];
h.Color = [0.5 0.5 0.5];
h.LineStyle = '-';
h.LineWidth = 2;
leg(2) = h;
hold on
% h = errorbarxy(avg_oric(oric_ind), HR_av(oric_ind), sem_oric(oric_ind), sem_oric(oric_ind), HR_av(oric_ind) - ci95_HR_av((oric_ind),1)', ci95_HR_av((oric_ind),2)' - HR_av(oric_ind), {'oc', 'c', 'c'});
h = plot(avg_oric(oric_ind), FR_oricpvv(oric_ind),'co');
h.MarkerFaceColor  = 'c';
h.LineStyle = '-';
h.LineWidth = 2;
leg(3) = h;
h = plot(avg_oric(oric_ind), FR_oricpaa(oric_ind),'co');
h.MarkerFaceColor  = [0 0.75 0.75];
h.Color = [0 0.75 0.75];
h.LineStyle = '-';
h.LineWidth = 2;
leg(4) = h;
x_tick = [0 x_axis_ori];
x_tick_label = [0 x_axis_ori];
y_tick = y_axis_hr;
figXAxis(h.Parent,'Orientation change (deg)',[0 100],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_tick);
figAxForm(h.Parent)
legend(leg,{'vis->vis->vis';'aud->aud->vis';'vis->vis->inv vis';'aud->aud->inv vis'},'location','northwest')

subplot 122
% h = errorbar(avg_trTime,HR_pv_trTime,HR_pv_trTime-ci_95_HR_pv_trTime(:,1)',ci_95_HR_pv_trTime(:,2)'-HR_pv_trTime,'ko');
h = plot(avg_trTime,HR_pv_trTime,'ko');
h.MarkerFaceColor = 'k';
hold on
% h = errorbar(avg_trTime,HR_pa_trTime,HR_pa_trTime-ci_95_HR_pa_trTime(:,1)',ci_95_HR_pa_trTime(:,2)'-HR_pa_trTime,'ko');
h = plot(avg_trTime,HR_pa_trTime,'ko');
h.MarkerFaceColor = [0.5 0.5 0.5];
hold on
% h = errorbar(avg_trTime,FR_pv_trTime,FR_pv_trTime-ci_95_FR_pv_trTime(:,1)',ci_95_FR_pv_trTime(:,2)'-FR_pv_trTime,'co');
h = plot(avg_trTime,FR_pv_trTime,'co');
h.MarkerFaceColor = [0 0.75 0.75];
h.Color = [0 0.75 0.75];
% h = errorbar(avg_trTime,FR_pa_trTime,FR_pa_trTime-ci_95_FR_pa_trTime(:,1)',ci_95_FR_pa_trTime(:,2)'-FR_pa_trTime,'co');
h = plot(avg_trTime,FR_pa_trTime,'co');
h.MarkerFaceColor = 'c';
hold on
x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(h.Parent)
title('visual trials')

print(fullfile(fnout,'all_expt_depPrevTrial'),'-dpdf','-fillpage')