    
% visInd_all = find(visTargets_all > 1);
% cVisInd_all = find(invVisTargets_all > 1);
% prevVis = find(prevTrial_all == 1);
% prevAud = find(prevTrial_all == 0);

nTimes = 6;
trialTime_edges = [0:max(trLengthMs_all)/nTimes:max(trLengthMs_all)]+100;

[h_time, ~, bin_time] = histcounts(trLengthMs_all,trialTime_edges);

time_bins = unique(bin_time);
avg_trTime = zeros(1,length(time_bins));
hits_v_trTime = zeros(1,length(time_bins));
misses_v_trTime = zeros(1,length(time_bins));
hits_av_trTime = zeros(1,length(time_bins));
misses_av_trTime = zeros(1,length(time_bins));
% hits_pv_trTime = zeros(1,length(time_bins));
% misses_pv_trTime = zeros(1,length(time_bins));
% hits_pa_trTime = zeros(1,length(time_bins));
% misses_pa_trTime = zeros(1,length(time_bins));

for ibin = 1:length(time_bins)
    ind_all = find(bin_time == ibin);
    avg_trTime(ibin) = mean(trLengthMs_all(ind_all),2);
    
    ind_v = bin_time == ibin & trType_all == visual;
    ind_av = bin_time == ibin & trType_all == auditory;
    
%     v_pv = intersect(v,prevVis);
%     v_pa = intersect(v,prevAud);
%     c_pv = intersect(c,prevVis);
%     c_pa = intersect(c,prevAud);


    hits_v_trTime(ibin) = sum(successIx_all(ind_v));
    misses_v_trTime(ibin) = sum(missedIx_all(ind_v));
    hits_av_trTime(ibin) = sum(successIx_all(ind_av));
    misses_av_trTime(ibin) = sum(missedIx_all(ind_av));
        
    ind_v = bin_time == ibin & trType_all == visual & targets_all < thresh;
    ind_av = bin_time == ibin & trType_all == auditory & targets_all < thresh;
    
    hits_v_thresh_trTime(ibin) = sum(successIx_all(ind_v));
    misses_v_thresh_trTime(ibin) = sum(missedIx_all(ind_v));
    hits_av_thresh_trTime(ibin) = sum(successIx_all(ind_av));
    misses_av_thresh_trTime(ibin) = sum(missedIx_all(ind_av));
        
%     hits_pv_trTime(ibin) = sum(successIx_all(v_pv));
%     misses_pv_trTime(ibin) = sum(missedIx_all(v_pv));
%     
%     hits_pa_trTime(ibin) = sum(successIx_all(v_pa));
%     misses_pa_trTime(ibin) = sum(missedIx_all(v_pa));
%     
%     invHit_trTime(ibin) = sum(invHitIx_all(c));
%     invMiss_trTime(ibin) = sum(invMissIx_all(c));
%     
%     invHit_pv_trTime(ibin) = sum(invHitIx_all(c_pv));
%     invMiss_pv_trTime(ibin) = sum(invMissIx_all(c_pv));
%     
%     invHit_pa_trTime(ibin) = sum(invHitIx_all(c_pa));
%     invMiss_pa_trTime(ibin) = sum(invMissIx_all(c_pa));
    
end
[HR_v_trTime ci_95_HR_v_trTime] = binofit(hits_v_trTime,hits_v_trTime+misses_v_trTime);
[HR_av_trTime ci_95_HR_av_trTime] = binofit(hits_av_trTime,hits_av_trTime+misses_av_trTime);

[HR_v_thresh_trTime ci_95_HR_v_thresh_trTime] = binofit(hits_v_thresh_trTime,hits_v_thresh_trTime+misses_v_thresh_trTime);
[HR_av_thresh_trTime ci_95_HR_av_thresh_trTime] = binofit(hits_av_thresh_trTime,hits_av_thresh_trTime+misses_av_thresh_trTime);

% [HR_pv_trTime ci_95_HR_pv_trTime] = binofit(hits_pv_trTime,hits_pv_trTime+misses_pv_trTime);
% [HR_pa_trTime ci_95_HR_pa_trTime] = binofit(hits_pa_trTime,hits_pa_trTime+misses_pa_trTime);
% 
% [FR_pv_trTime ci_95_FR_pv_trTime] = binofit(invHit_pv_trTime,invHit_pv_trTime+invMiss_pv_trTime);
% [FR_pa_trTime ci_95_FR_pa_trTime] = binofit(invHit_pa_trTime,invHit_pa_trTime+invMiss_pa_trTime);

xTick_ind = find(~isnan(avg_trTime));
