if hr_calc == 1
[h_all, ~, bin_all] = histcounts(targets, ori_edges);
[h_v, ~, bin_v] = histcounts(vTargets, ori_edges);
[h_av, ~, bin_av] = histcounts(avTargets, ori_edges);

ori_bins = unique(bin_v);
avg_ori = zeros(1,max(ori_bins,[],2));
sem_ori = zeros(1,max(ori_bins,[],2));
Hits_v = zeros(1,max(ori_bins,[],2));
Misses_v = zeros(1,max(ori_bins,[],2));
Hit_av = zeros(1,max(ori_bins,[],2));
Misses_av = zeros(1,max(ori_bins,[],2));
n_ori = zeros(2,max(ori_bins,[],2));
for ibin = 1:max(ori_bins,[],2)
    ind_all = find(bin_all == ibin);
    avg_ori(ibin) = mean(targets(:,ind_all),2);
    sem_ori(ibin) = std(targets(:,ind_all),[],2)./sqrt(length(ind_all));
    
    ind = bin_all == ibin & trType == visual;
    if sum(ind)>=minTrN_expt
        Hits_v(ibin) = sum(successIx(ind),2);
        Misses_v(ibin) = sum(missedIx(ind),2);
    end
    ind = bin_all == ibin & trType == auditory;
    if sum(ind)>=minTrN_expt
        Hit_av(ibin) = sum(successIx(ind),2);
        Misses_av(ibin) = sum(missedIx(ind),2);
    end
end
[HR_v, ci95_HR_v] = binofit(Hits_v, Misses_v + Hits_v);
[HR_av, ci95_HR_av] = binofit(Hit_av, Misses_av + Hit_av);
n_ori = [Hits_v+Misses_v; Hit_av+Misses_av];

%%
elseif hr_calc == 2 % across mice
[h_all, ~, bin_all] = histcounts(targets_all, ori_edges);
[h_v, ~, bin_v] = histcounts(vTargets_all, ori_edges);
[h_av, ~, bin_av] = histcounts(avTargets_all, ori_edges);
    

[h_all, ~, bin_all] = histcounts(targets_all, ori_edges); 
% [h_v, ~, bin_v] = histcounts(vTargets_all, ori_edges); 
% [h_av, ~, bin_av] = histcounts(avTargets_all, ori_edges);

ori_bins = unique(bin_v);
avg_ori = zeros(1,max(ori_bins,[],2));
sem_ori = zeros(1,max(ori_bins,[],2));
Hits_v = zeros(1,max(ori_bins,[],2));
Misses_v = zeros(1,max(ori_bins,[],2));
Hit_av = zeros(1,max(ori_bins,[],2));
Misses_av = zeros(1,max(ori_bins,[],2));

% %prev trial
% Hits_vpv = zeros(1,max(ori_bins,[],2));
% Misses_vpv = zeros(1,max(ori_bins,[],2));
% Hits_vpa = zeros(1,max(ori_bins,[],2));
% Misses_vpa = zeros(1,max(ori_bins,[],2));
% Hits_avpv = zeros(1,max(ori_bins,[],2));
% Misses_avpv = zeros(1,max(ori_bins,[],2));
% Hits_avpa = zeros(1,max(ori_bins,[],2));
% Misses_avpa = zeros(1,max(ori_bins,[],2));
% %prev 2 trials
% Hits_vpvv = zeros(1,max(ori_bins,[],2));
% Misses_vpvv = zeros(1,max(ori_bins,[],2));
% Hits_vpaa = zeros(1,max(ori_bins,[],2));
% Misses_vpaa = zeros(1,max(ori_bins,[],2));
% Hits_avpvv = zeros(1,max(ori_bins,[],2));
% Misses_avpvv = zeros(1,max(ori_bins,[],2));
% Hits_avpaa = zeros(1,max(ori_bins,[],2));
% Misses_avpaa = zeros(1,max(ori_bins,[],2));

n_ori = zeros(2,max(ori_bins,[],2));
for ibin = 1:max(ori_bins,[],2)
    ind_all = find(bin_all == ibin);
    avg_ori(ibin) = mean(targets_all(ind_all),2);
    sem_ori(ibin) = std(targets_all(:,ind_all),[],2)./sqrt(length(ind_all));
    
    ind = bin_all == ibin & trType_all == visual;
    
%     indpv = find(bin_v == ibin & prevTrial_all == 1) ;
%     indpa = find(bin_v == ibin & prevTrial_all == 0) ;
%     indpvv = find(bin_v == ibin & prevTrial_all == 1 & prev2Trial_all == 1) ;
%     indpaa = find(bin_v == ibin & prevTrial_all == 0 & prev2Trial_all == 0) ;
    if sum(ind)>=minTrN_all
        Hits_v(ibin) = sum(successIx_all(ind),2);
        Misses_v(ibin) = sum(missedIx_all(ind),2);
        
%         Hits_vpv(ibin) = sum(successIx_all(indpv),2);
%         Misses_vpv(ibin) = sum(missedIx_all(indpv),2);
%         Hits_vpa(ibin) = sum(successIx_all(indpa),2);
%         Misses_vpa(ibin) = sum(missedIx_all(indpa),2);
        
%         Hits_vpvv(ibin) = sum(successIx_all(indpvv),2);
%         Misses_vpvv(ibin) = sum(missedIx_all(indpvv),2);
%         Hits_vpaa(ibin) = sum(successIx_all(indpaa),2);
%         Misses_vpaa(ibin) = sum(missedIx_all(indpaa),2);
    end
    ind = bin_all == ibin & trType_all == auditory;
%     indcpv = find(bin_av == ibin & prevTrial_all == 1) ;
%     indcpa = find(bin_av == ibin & prevTrial_all == 0) ;
%     indcpvv = find(bin_av == ibin & prevTrial_all == 1 & prev2Trial_all == 1) ;
%     indcpaa = find(bin_av == ibin & prevTrial_all == 0 & prev2Trial_all == 0) ;
    if sum(ind)>=minTrN_all
        Hit_av(ibin) = sum(successIx_all(ind),2);
        Misses_av(ibin) = sum(missedIx_all(ind),2);
        
%         Hits_avpv(ibin) = sum(invHitIx_all(indcpv),2);
%         Misses_avpv(ibin) = sum(invMissIx_all(indcpv),2);
%         Hits_avpa(ibin) = sum(invHitIx_all(indcpa),2);
%         Misses_avpa(ibin) = sum(invMissIx_all(indcpa),2);        
        
%         Hits_avpvv(ibin) = sum(invHitIx_all(indcpvv),2);
%         Misses_avpvv(ibin) = sum(invMissIx_all(indcpvv),2);
%         Hits_avpaa(ibin) = sum(invHitIx_all(indcpaa),2);
%         Misses_avpaa(ibin) = sum(invMissIx_all(indcpaa),2);
        
    end
end
[HR_v, ci95_HR_v] = binofit(Hits_v, Misses_v + Hits_v);
[HR_av, ci95_HR_av] = binofit(Hit_av, Misses_av + Hit_av);

% [HR_vpv, ci95_HR_vpv] = binofit(Hits_vpv, Misses_vpv + Hits_vpv);
% [HR_vpa, ci95_HR_vpa] = binofit(Hits_vpa, Misses_vpa + Hits_vpa);
% 
% [HR_vpvv, ci95_HR_vpvv] = binofit(Hits_vpvv, Misses_vpvv + Hits_vpvv);
% [HR_vpaa, ci95_HR_vpaa] = binofit(Hits_vpaa, Misses_vpaa + Hits_vpaa);
% 
% [HR_avpv, ci95_HR_avpv] = binofit(Hits_avpv, Misses_avpv + Hits_avpv);
% [HR_avpa, ci95_HR_avpa] = binofit(Hits_avpa, Misses_avpa + Hits_avpa);
% 
% [HR_avpvv, ci95_HR_avpvv] = binofit(Hits_avpvv, Misses_avpvv + Hits_avpvv);
% [HR_avpaa, ci95_HR_avpaa] = binofit(Hits_avpaa, Misses_avpaa + Hits_avpaa);

n_ori = [Hits_v+Misses_v; Hit_av+Misses_av];
    
end