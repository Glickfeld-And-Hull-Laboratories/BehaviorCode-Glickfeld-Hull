%% run beginning of plotAVCatchsummary.m first.
%% what are the react times for success trials?
RTs = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        RTs = [RTs mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(hits).reactTimes];
    end
end

TCs_val_inv = cat(2,mean(tc_val_all,2)-mean(mean(tc_val_all(pre_win,:),2)),mean(tc_inv_all,2)-mean(mean(tc_inv_all(pre_win,:),2)));

[nRTs binRTs] = hist(RTs,10);

figure;
[ax,h_rt,h_tc] = plotyy(binRTs,nRTs,ttMs,TCs_val_inv,'bar','plot');
ax(2).YLim = [-0.0025 0.03];
ax(1).YLabel.String = 'n trials';
ax(2).YLabel.String = 'dF/F';
ax(1).XLim = [-10 20]/(cycTime/cycTimeMs);
ax(2).XLim = [-10 20]/(cycTime/cycTimeMs);
h_rt.FaceColor = [0.75 0.75 0.75];
h_tc(1).Color = 'k';
h_tc(2).Color = 'c';
hold on
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title([num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']);


print([fnout 'tc_val_inv_withRT.pdf'], '-dpdf')

%% react times for vis, aud, and catch successes

RT_vis_h = [];
resp_vis_h_tc = mean(tc_hvsfa_all,2)-mean(mean(tc_hvsfa_all(pre_event_frames,:),2));
RT_aud_h = [];
resp_aud_h_tc_all = [];
RT_inv_h = [];
resp_inv_h_tc = mean(tc_favsh_all,2)-mean(mean(tc_favsh_all(pre_event_frames,:),2));

resp_vis_m_tc = mean(tc_mvscr_all,2)-mean(mean(tc_mvscr_all(pre_event_frames,:),2));
resp_aud_m_tc_all = [];
resp_inv_m_tc = mean(tc_crvsm_all,2)-mean(mean(tc_crvsm_all(pre_event_frames,:),2));

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if mouse(imouse).expt(iexp).info.isCatch
        RT_vis_h = [RT_vis_h mouse(imouse).expt(iexp).align(2).av(1).outcome(1).reactTimes];
        RT_aud_h = [RT_aud_h mouse(imouse).expt(iexp).align(2).av(2).outcome(1).reactTimes];
        RT_inv_h = [RT_inv_h mouse(imouse).expt(iexp).align(3).av(1).FAreacttime];
        
        resp_aud_h_tc_all = cat(2,resp_aud_h_tc_all,mean(mouse(imouse).expt(iexp).align(2).av(2).outcome(1).resp,3));
        resp_aud_m_tc_all = cat(2,resp_aud_m_tc_all,mean(mouse(imouse).expt(iexp).align(2).av(2).outcome(2).resp,3));
        end
    end
end
resp_aud_h_tc = mean(resp_aud_h_tc_all,2)-mean(mean(resp_aud_h_tc_all(pre_win,:),2));
resp_aud_m_tc = nanmean(resp_aud_m_tc_all,2)-nanmean(nanmean(resp_aud_m_tc_all(pre_win,:),2));

[nRT_vis binRT_vis] = hist(RT_vis_h,10);
[nRT_aud binRT_aud] = hist(RT_aud_h,10);
[nRT_inv binRT_inv] = hist(RT_inv_h,10);


% plot vis RTs and timecourse
%vis trials
RTvsRespFig = figure;
subplot(3,2,1)
[ax,h_tc,h_rt] = plotyy(ttMs,resp_vis_h_tc,binRT_vis,nRT_vis,'plot','bar');
ax(1).YLim = [-0.0025 0.03];
ax(2).YLabel.String = 'n trials';
ax(1).YLabel.String = 'dF/F';
ax(1).YTick = [0 0.01 0.02 0.03];
ax(2).XLim = [-10 20]/(cycTime/cycTimeMs);
ax(1).XLim = [-10 20]/(cycTime/cycTimeMs);
h_rt.FaceColor = [0.75 0.75 0.75];
h_tc.Color = 'g';
hold on
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['visual trials - all valid hits, ' num2str(size(tc_hvsfa_all,2)) ' cells']);

subplot(3,2,2)
plot(ttMs,resp_vis_m_tc,'r')
hold on
xlim([-10 20]/(cycTime/cycTimeMs));
ylim([-0.0025 0.03]);
ax = gca;
ax.YTick = [0 0.01 0.02 0.03];
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['visual trials - miss match to CR, ' num2str(size(tc_mvscr_all,2)) ' cells']);


%invalid vis trials
figure(RTvsRespFig);
subplot(3,2,3)
[ax,h_tc,h_rt] = plotyy(ttMs,resp_inv_h_tc,binRT_inv,nRT_inv,'plot','bar');
ax(1).YLim = [-0.0025 0.03];
ax(2).YLabel.String = 'n trials';
ax(1).YLabel.String = 'dF/F';
ax(1).YTick = [0 0.01 0.02 0.03];
ax(2).XLim = [-10 20]/(cycTime/cycTimeMs);
ax(1).XLim = [-10 20]/(cycTime/cycTimeMs);
h_rt.FaceColor = [0.75 0.75 0.75];
h_tc.Color = 'c';
hold on
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['visual trials - all invalid hits, ' num2str(size(resp_inv_all,2)) ' cells']);

subplot(3,2,4)
plot(ttMs,resp_inv_m_tc,'b')
hold on
xlim([-10 20]/(cycTime/cycTimeMs));
ylim([-0.0025 0.03]);
ax = gca;
ax.YTick = [0 0.01 0.02 0.03];
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['visual trials - CR match to m, ' num2str(size(tc_crvsm_all,2)) ' cells']);

%aud trials
figure(RTvsRespFig);
subplot(3,2,5)
[ax,h_tc,h_rt] = plotyy(ttMs,resp_aud_h_tc,binRT_aud,nRT_aud,'plot','bar');
ax(1).YLim = [-0.0025 0.03];
ax(2).YLabel.String = 'n trials';
ax(1).YLabel.String = 'dF/F';
ax(1).YTick = [0 0.01 0.02 0.03];
ax(2).XLim = [-10 20]/(cycTime/cycTimeMs);
ax(1).XLim = [-10 20]/(cycTime/cycTimeMs);
h_rt.FaceColor = [0.75 0.75 0.75];
h_tc.Color = 'k';
hold on
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['auditory trials - all hits, ' num2str(size(resp_aud_h_tc_all,2)) ' cells']);

subplot(3,2,6)
plot(ttMs,resp_aud_m_tc,'m')
hold on
xlim([-10 20]/(cycTime/cycTimeMs));
ylim([-0.0025 0.03]);
ax = gca;
ax.YTick = [0 0.01 0.02 0.03];
vline(0,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title(['aud trials - all miss, ' num2str(size(resp_aud_m_tc_all,2)) ' cells']);


print([fnout 'tc_A-V-inv_withRT.pdf'], '-dpdf')