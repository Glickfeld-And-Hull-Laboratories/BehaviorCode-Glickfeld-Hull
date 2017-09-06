clear all
close all
rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename_audCtrl, [], rc.indexTextCols);
datasetStr = ['awFSAVdatasets_audControl'];
fnout = fullfile(rc.caOutputDir,datasetStr,'behavior');
if ~exist(fnout,'dir')
    mkdir(fnout)
end
%%
minTrN_expt = 10;
minTrN_all = 10;
thresh = 30;

nexp = xd.nRows;
mice = unique(xd.Subject);
nmice = length(mice);
colors_mice = parula(nmice);
visual = 0;
auditory = 1;

doExptStruct = 1;
if doExptStruct
    bxExptStruct_FSaudcon_train
else
    load(fullfile(fnout,'bxExpMat_training.mat'));
end

amp_edges = [0.003 0.007 0.02 0.05 0.14 0.37 1];
% amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

%%
ms = struct;
for im = 1:nmice
    ms(im).sn = mice(im);
end
targets_all = [];
vTargets_all = [];
avTargets_all = [];
successIx_all = [];
missedIx_all = [];
trLengthMs_all = [];
trType_all = [];
prevTrial_all = [];
prev2Trial_all = [];

vHR_expt = zeros(1,nexp);
avHR_expt = zeros(1,nexp);

doPlot = 0;
%% analyze behavior data from each experiment
% hr_calc = 1;
for iexp = 1:nexp
%     title_str = [bxExp(iexp).sn '-' expt(iexp).date];
    sn = bxExp(iexp).sn;
    sn_ind = find(mice == sn);
    cycLengthMs = double(bxExp(iexp).tOn + bxExp(iexp).tOff);
    
    targets = chop(double(bxExp(iexp).allTargets),2);
    trType = bxExp(iexp).trType;
    
    trType_shift = [NaN trType];
    prevTrial = trType_shift(1:length(trType));
    trType_shift = [NaN NaN trType];
    prev2Trial = trType_shift(1:length(trType));
   
    successIx = double(bxExp(iexp).sIx);
    missedIx = double(bxExp(iexp).mIx);
    nCyc = double(bxExp(iexp).trLength);
    trLengthMs = nCyc.*cycLengthMs;
    
    vTargets = targets(trType == visual);
    avTargets = targets(trType == auditory);
    
    if any(strcmp(fieldnames(ms(sn_ind)),'targets_all'))
        ms(sn_ind).targets_all = cat(2,ms(sn_ind).targets_all,targets);
        ms(sn_ind).vTargets_all = cat(2,ms(sn_ind).vTargets_all,vTargets);
        ms(sn_ind).avTargets_all = cat(2,ms(sn_ind).avTargets_all,avTargets);
        ms(sn_ind).successIx_all = cat(2,ms(sn_ind).successIx_all,successIx);
        ms(sn_ind).missedIx_all = cat(2,ms(sn_ind).missedIx_all,missedIx);
        ms(sn_ind).trLengthMs_all = cat(2,ms(sn_ind).trLengthMs_all,trLengthMs);
        ms(sn_ind).trType_all = cat(2,ms(sn_ind).trType_all,trType);
        ms(sn_ind).prevTrial_all = cat(2,ms(sn_ind).prevTrial_all,prevTrial);
        ms(sn_ind).prev2Trial_all = cat(2,ms(sn_ind).prev2Trial_all,prev2Trial);
    else
        ms(sn_ind).targets_all = targets;
        ms(sn_ind).vTargets_all = vTargets;
        ms(sn_ind).avTargets_all = avTargets;
        ms(sn_ind).successIx_all = successIx;
        ms(sn_ind).missedIx_all = missedIx;
        ms(sn_ind).trLengthMs_all = trLengthMs;
        ms(sn_ind).trType_all = trType;
        ms(sn_ind).prevTrial_all = prevTrial;
        ms(sn_ind).prev2Trial_all = prev2Trial;        
    end
    targets_all = cat(2,targets_all,targets);
    vTargets_all = cat(2,vTargets_all,vTargets);
    avTargets_all = cat(2,avTargets_all,avTargets);
    successIx_all = cat(2,successIx_all,successIx);
    missedIx_all = cat(2,missedIx_all,missedIx);
    trLengthMs_all = cat(2,trLengthMs_all,trLengthMs);
    trType_all = cat(2,trType_all,trType);
    prevTrial_all = cat(2,prevTrial_all,prevTrial);
    prev2Trial_all = cat(2,prev2Trial_all,prev2Trial);
    
%     bxSumExptHR_FSaudcon_img
    
    vHR_expt(iexp) = sum(successIx & trType == visual)./(sum(successIx & trType == visual) + sum(missedIx & trType == visual));
    avHR_expt(iexp) = sum(successIx & trType == auditory)./(sum(successIx & trType == auditory) + sum(missedIx & trType == auditory));
end

%% analyze behavior data across all experiments
hr_calc = 2;

bxSumExptHR_FSaudcon_img

v_ind = ~isnan(HR_v);
av_ind = ~isnan(HR_av);
x_axis_ori = fliplr(chop(spoc(90,1,4),2));

% hit rate - valid and invalid visual trials
y_axis_hr = 0:0.2:1;
all_fig = figure; setFigParams4Print('portrait')
suptitle('all training expt')
figure(all_fig);
sp2 = subplot(2,2,1);
h = errorbarxy(avg_ori(v_ind), HR_v(v_ind), sem_ori(v_ind), sem_ori(v_ind), HR_v(v_ind) - ci95_HR_v((v_ind),1)', ci95_HR_v((v_ind),2)' - HR_v(v_ind), {'ok', 'k', 'k'});
h.hMain.MarkerFaceColor  = 'k';
h.hMain.LineStyle = '-';
h.hMain.LineWidth = 3;
hold on
h = errorbarxy(avg_ori(av_ind), HR_av(av_ind), sem_ori(av_ind), sem_ori(av_ind), HR_av(av_ind) - ci95_HR_av((av_ind),1)', ci95_HR_av((av_ind),2)' - HR_av(av_ind), {'oc', 'c', 'c'});
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
title({['vis = ' num2str(sum(n_ori(1,:),2))], ['vis+aud = ' num2str(sum(n_ori(2,:),2))]});
% axis square

% hit rate by trial length
bxSumTiming_FSaudcon_img

v_ind = ~isnan(HR_v_trTime);
av_ind = ~isnan(HR_av_trTime);

figure(all_fig)
sp3 = subplot(2,2,2);
h = errorbar(avg_trTime(v_ind),HR_v_trTime(v_ind),HR_v_trTime(v_ind)-ci_95_HR_v_trTime((v_ind),1)',ci_95_HR_v_trTime((v_ind),2)'-HR_v_trTime(v_ind),'ko');
h.MarkerFaceColor = 'k';
h.LineStyle = '-';
h.LineWidth = 1;
leg(1) = h;
hold on
h = errorbar(avg_trTime(av_ind),HR_av_trTime(av_ind),HR_av_trTime(av_ind)-ci_95_HR_av_trTime((av_ind),1)',ci_95_HR_av_trTime((av_ind),2)'-HR_av_trTime(av_ind),'co');
h.MarkerFaceColor = 'c';
h.LineStyle = '-';
h.LineWidth = 1;
leg(2) = h;
hold on

x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(sp3)
title('binned by trial length')
l = legend(leg,{'vis';'vis+aud'},'location','southeastoutside');
title(l,'trial type')

% HR x trial length x trial difficulty
figure(all_fig)
sp3 = subplot(2,2,3);
h = errorbar(avg_trTime(v_ind),HR_v_thresh_trTime(v_ind),HR_v_thresh_trTime(v_ind)-ci_95_HR_v_thresh_trTime((v_ind),1)',ci_95_HR_v_thresh_trTime((v_ind),2)'-HR_v_thresh_trTime(v_ind),'ko');
h.MarkerFaceColor = 'k';
h.LineStyle = '-';
h.LineWidth = 1;
leg(1) = h;
hold on
h = errorbar(avg_trTime(av_ind),HR_av_thresh_trTime(av_ind),HR_av_thresh_trTime(av_ind)-ci_95_HR_av_thresh_trTime((av_ind),1)',ci_95_HR_av_thresh_trTime((av_ind),2)'-HR_av_thresh_trTime(av_ind),'co');
h.MarkerFaceColor = 'c';
h.LineStyle = '-';
h.LineWidth = 1;
leg(2) = h;
hold on

x_tick = [0 chop(avg_trTime,2)];
x_tick_label = [0 chop(avg_trTime./1000,2)];
figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
figAxForm(sp3)
title({'binned by trial length';['target < ' num2str(thresh)]})
l = legend(leg,{'vis';'vis+aud'},'location','southeastoutside');
title(l,'trial type')

figure(all_fig)
print(fullfile(fnout,'training_all_expt'),'-dpdf','-fillpage')
%%
% % %% behavior dep on previous trial
% % % hit rate - vis trials sorted by preceeding trial type
% % 
% % figure; setFigParams4Print('landscape')
% % subplot 121
% % % h = errorbarxy(avg_ori(ori_ind), HR_oripv(ori_ind), sem_ori(ori_ind), sem_ori(ori_ind), HR_oripv(ori_ind) - ci95_HR_oripv((ori_ind),1)', ci95_HR_oripv((ori_ind),2)' - HR_oripv(ori_ind), {'ok', 'k', 'k'});
% % h = plot(avg_ori(v_ind), HR_oripvv(v_ind),'ko');
% % h.MarkerFaceColor  = 'k';
% % h.LineStyle = '-';
% % h.LineWidth = 2;
% % leg(1) = h;
% % hold on
% % % h = errorbarxy(avg_ori(ori_ind), HR_oripa(ori_ind), sem_ori(ori_ind), sem_ori(ori_ind), HR_oripa(ori_ind) - ci95_HR_oripa((ori_ind),1)', ci95_HR_oripa((ori_ind),2)' - HR_oripa(ori_ind), {'ok', 'k', 'k'});
% % h = plot(avg_ori(v_ind), HR_oripaa(v_ind),'ko');
% % h.MarkerFaceColor  = [0.5 0.5 0.5];
% % h.Color = [0.5 0.5 0.5];
% % h.LineStyle = '-';
% % h.LineWidth = 2;
% % leg(2) = h;
% % hold on
% % % h = errorbarxy(avg_ori(oric_ind), HR_av(oric_ind), sem_ori(oric_ind), sem_ori(oric_ind), HR_av(oric_ind) - ci95_HR_av((oric_ind),1)', ci95_HR_av((oric_ind),2)' - HR_av(oric_ind), {'oc', 'c', 'c'});
% % h = plot(avg_ori(av_ind), FR_oricpvv(av_ind),'co');
% % h.MarkerFaceColor  = 'c';
% % h.LineStyle = '-';
% % h.LineWidth = 2;
% % leg(3) = h;
% % h = plot(avg_ori(av_ind), FR_oricpaa(av_ind),'co');
% % h.MarkerFaceColor  = [0 0.75 0.75];
% % h.Color = [0 0.75 0.75];
% % h.LineStyle = '-';
% % h.LineWidth = 2;
% % leg(4) = h;
% % x_tick = [0 x_axis_ori];
% % x_tick_label = [0 x_axis_ori];
% % y_tick = y_axis_hr;
% % figXAxis(h.Parent,'Orientation change (deg)',[0 100],x_tick,x_tick_label);
% % figYAxis(h.Parent,'HR',[0 1.1],y_tick);
% % figAxForm(h.Parent)
% % legend(leg,{'vis->vis->vis';'aud->aud->vis';'vis->vis->inv vis';'aud->aud->inv vis'},'location','northwest')
% % 
% % subplot 122
% % % h = errorbar(avg_trTime,HR_pv_trTime,HR_pv_trTime-ci_95_HR_pv_trTime(:,1)',ci_95_HR_pv_trTime(:,2)'-HR_pv_trTime,'ko');
% % h = plot(avg_trTime,HR_pv_trTime,'ko');
% % h.MarkerFaceColor = 'k';
% % hold on
% % % h = errorbar(avg_trTime,HR_pa_trTime,HR_pa_trTime-ci_95_HR_pa_trTime(:,1)',ci_95_HR_pa_trTime(:,2)'-HR_pa_trTime,'ko');
% % h = plot(avg_trTime,HR_pa_trTime,'ko');
% % h.MarkerFaceColor = [0.5 0.5 0.5];
% % hold on
% % % h = errorbar(avg_trTime,FR_pv_trTime,FR_pv_trTime-ci_95_FR_pv_trTime(:,1)',ci_95_FR_pv_trTime(:,2)'-FR_pv_trTime,'co');
% % h = plot(avg_trTime,FR_pv_trTime,'co');
% % h.MarkerFaceColor = [0 0.75 0.75];
% % h.Color = [0 0.75 0.75];
% % % h = errorbar(avg_trTime,FR_pa_trTime,FR_pa_trTime-ci_95_FR_pa_trTime(:,1)',ci_95_FR_pa_trTime(:,2)'-FR_pa_trTime,'co');
% % h = plot(avg_trTime,FR_pa_trTime,'co');
% % h.MarkerFaceColor = 'c';
% % hold on
% % x_tick = [0 chop(avg_trTime,2)];
% % x_tick_label = [0 chop(avg_trTime./1000,2)];
% % figXAxis(h.Parent,'trial length (s)',[0 trialTime_edges(end)+500],x_tick,x_tick_label);
% % figYAxis(h.Parent,'HR',[0 1.1],y_axis_hr)
% % figAxForm(h.Parent)
% % title('visual trials')
% % 
% % print(fullfile(fnout,'all_expt_depPrevTrial'),'-dpdf','-fillpage')