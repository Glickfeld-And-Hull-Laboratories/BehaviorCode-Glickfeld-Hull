    
visInd_all = find(visTargets_all > 1);
cVisInd_all = find(invVisTargets_all > 1);
prevVis = find(prevTrial_all == 1);
prevAud = find(prevTrial_all == 0);

nTimes = 6;
trialTime_edges = [0:max(trLengthMs_all)/nTimes:max(trLengthMs_all)]+100;

[h_time, edge, bin_time] = histcounts(trLengthMs_all,trialTime_edges);

time_bins = unique(bin_time);
avg_trTime = zeros(1,length(time_bins));
hits_trTime = zeros(1,length(time_bins));
misses_trTime = zeros(1,length(time_bins));
hits_pv_trTime = zeros(1,length(time_bins));
misses_pv_trTime = zeros(1,length(time_bins));
hits_pa_trTime = zeros(1,length(time_bins));
misses_pa_trTime = zeros(1,length(time_bins));
invHit_trTime = zeros(1,length(time_bins));
invMiss_trTime = zeros(1,length(time_bins));
invHit_pv_trTime = zeros(1,length(time_bins));
invMiss_pv_trTime = zeros(1,length(time_bins));
invHit_pa_trTime = zeros(1,length(time_bins));
invMiss_pa_trTime = zeros(1,length(time_bins));

for ibin = 1:length(time_bins)
    ind_v = find(bin_time == time_bins(ibin));
    v = intersect(visInd_all,ind_v);    
    c = intersect(cVisInd_all,ind_v); 
    
    v_pv = intersect(v,prevVis);
    v_pa = intersect(v,prevAud);
    c_pv = intersect(c,prevVis);
    c_pa = intersect(c,prevAud);

    avg_trTime(ibin) = mean(trLengthMs_all(ind_v));

    hits_trTime(ibin) = sum(successIx_all(v));
    misses_trTime(ibin) = sum(missedIx_all(v));
        
    hits_pv_trTime(ibin) = sum(successIx_all(v_pv));
    misses_pv_trTime(ibin) = sum(missedIx_all(v_pv));
    
    hits_pa_trTime(ibin) = sum(successIx_all(v_pa));
    misses_pa_trTime(ibin) = sum(missedIx_all(v_pa));
    
    invHit_trTime(ibin) = sum(invHitIx_all(c));
    invMiss_trTime(ibin) = sum(invMissIx_all(c));
    
    invHit_pv_trTime(ibin) = sum(invHitIx_all(c_pv));
    invMiss_pv_trTime(ibin) = sum(invMissIx_all(c_pv));
    
    invHit_pa_trTime(ibin) = sum(invHitIx_all(c_pa));
    invMiss_pa_trTime(ibin) = sum(invMissIx_all(c_pa));
    
end
[HR_trTime ci_95_HR_trTime] = binofit(hits_trTime,hits_trTime+misses_trTime);
[FR_trTime ci_95_FR_trTime] = binofit(invHit_trTime,invHit_trTime+invMiss_trTime);

[HR_pv_trTime ci_95_HR_pv_trTime] = binofit(hits_pv_trTime,hits_pv_trTime+misses_pv_trTime);
[HR_pa_trTime ci_95_HR_pa_trTime] = binofit(hits_pa_trTime,hits_pa_trTime+misses_pa_trTime);

[FR_pv_trTime ci_95_FR_pv_trTime] = binofit(invHit_pv_trTime,invHit_pv_trTime+invMiss_pv_trTime);
[FR_pa_trTime ci_95_FR_pa_trTime] = binofit(invHit_pa_trTime,invHit_pa_trTime+invMiss_pa_trTime);

xTick_ind = find(~isnan(avg_trTime));
