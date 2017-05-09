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
    if length(ind)>=minTrN_expt
        Hits_v(ibin) = sum(successIx(ind),2);
        Misses_v(ibin) = sum(missedIx(ind),2);
    end
    ind = bin_all == ibin & trType == auditory;
    if length(ind)>=minTrN_expt
        Hit_av(ibin) = sum(successIx(ind),2);
        Misses_av(ibin) = sum(missedIx(ind),2);
    end
end
[HR_v, ci95_HR_v] = binofit(Hits_v, Misses_v + Hits_v);
[HR_av, ci95_HR_av] = binofit(Hit_av, Misses_av + Hit_av);
n_ori = [Hits_v+Misses_v; Hit_av+Misses_av];

%%
elseif hr_calc == 2 % across mice
    
  
[h_v, edge, bin_v] = histcounts(visTargets_all, ori_edges);
[h_av, edge, bin_av] = histcounts(invVisTargets_all, ori_edges);

ori_bins = unique(bin_v);
avg_ori = zeros(1,max(ori_bins,[],2));
sem_ori = zeros(1,max(ori_bins,[],2));
Hits_v = zeros(1,max(ori_bins,[],2));
Misses_v = zeros(1,max(ori_bins,[],2));

%prev trial
Hits_oripv = zeros(1,max(ori_bins,[],2));
Misses_oripv = zeros(1,max(ori_bins,[],2));
Hits_oripa = zeros(1,max(ori_bins,[],2));
Misses_oripa = zeros(1,max(ori_bins,[],2));
Hits_oricpv = zeros(1,max(ori_bins,[],2));
Misses_oricpv = zeros(1,max(ori_bins,[],2));
Hits_oricpa = zeros(1,max(ori_bins,[],2));
Misses_oricpa = zeros(1,max(ori_bins,[],2));
%prev 2 trials
Hits_oripvv = zeros(1,max(ori_bins,[],2));
Misses_oripvv = zeros(1,max(ori_bins,[],2));
Hits_oripaa = zeros(1,max(ori_bins,[],2));
Misses_oripaa = zeros(1,max(ori_bins,[],2));
Hits_oricpvv = zeros(1,max(ori_bins,[],2));
Misses_oricpvv = zeros(1,max(ori_bins,[],2));
Hits_oricpaa = zeros(1,max(ori_bins,[],2));
Misses_oricpaa = zeros(1,max(ori_bins,[],2));

avg_oric = zeros(1,max(ori_bins,[],2));
sem_oric = zeros(1,max(ori_bins,[],2));
Hit_av = zeros(1,max(ori_bins,[],2));
Misses_av = zeros(1,max(ori_bins,[],2));
n_ori = zeros(2,max(ori_bins,[],2));
for ibin = 1:max(ori_bins,[],2)
    ind = find(bin_v == ibin);
    indpv = find(bin_v == ibin & prevTrial_all == 1) ;
    indpa = find(bin_v == ibin & prevTrial_all == 0) ;
    indpvv = find(bin_v == ibin & prevTrial_all == 1 & prev2Trial_all == 1) ;
    indpaa = find(bin_v == ibin & prevTrial_all == 0 & prev2Trial_all == 0) ;
    if length(ind)>=10
        avg_ori(ibin) = mean(visTargets_all(:,ind),2);
        sem_ori(ibin) = std(visTargets_all(:,ind),[],2)./sqrt(length(ind));
        Hits_v(ibin) = sum(successIx_all(ind),2);
        Misses_v(ibin) = sum(missedIx_all(ind),2);
        
        Hits_oripv(ibin) = sum(successIx_all(indpv),2);
        Misses_oripv(ibin) = sum(missedIx_all(indpv),2);
        Hits_oripa(ibin) = sum(successIx_all(indpa),2);
        Misses_oripa(ibin) = sum(missedIx_all(indpa),2);
        
        Hits_oripvv(ibin) = sum(successIx_all(indpvv),2);
        Misses_oripvv(ibin) = sum(missedIx_all(indpvv),2);
        Hits_oripaa(ibin) = sum(successIx_all(indpaa),2);
        Misses_oripaa(ibin) = sum(missedIx_all(indpaa),2);
    end
    indc = find(bin_av == ibin);
    indcpv = find(bin_av == ibin & prevTrial_all == 1) ;
    indcpa = find(bin_av == ibin & prevTrial_all == 0) ;
    indcpvv = find(bin_av == ibin & prevTrial_all == 1 & prev2Trial_all == 1) ;
    indcpaa = find(bin_av == ibin & prevTrial_all == 0 & prev2Trial_all == 0) ;
    if length(indc)>=10
        avg_oric(ibin) = mean(invVisTargets_all(:,indc),2);
        sem_oric(ibin) = std(invVisTargets_all(:,indc),[],2)./sqrt(length(indc));
        Hit_av(ibin) = sum(invHitIx_all(indc),2);
        Misses_av(ibin) = sum(invMissIx_all(indc),2);
        
        Hits_oricpv(ibin) = sum(invHitIx_all(indcpv),2);
        Misses_oricpv(ibin) = sum(invMissIx_all(indcpv),2);
        Hits_oricpa(ibin) = sum(invHitIx_all(indcpa),2);
        Misses_oricpa(ibin) = sum(invMissIx_all(indcpa),2);        
        
        Hits_oricpvv(ibin) = sum(invHitIx_all(indcpvv),2);
        Misses_oricpvv(ibin) = sum(invMissIx_all(indcpvv),2);
        Hits_oricpaa(ibin) = sum(invHitIx_all(indcpaa),2);
        Misses_oricpaa(ibin) = sum(invMissIx_all(indcpaa),2);
        
    end
end
[HR_v, ci95_HR_v] = binofit(Hits_v, Misses_v + Hits_v);
[HR_av, ci95_HR_av] = binofit(Hit_av, Misses_av + Hit_av);

[HR_oripv, ci95_HR_oripv] = binofit(Hits_oripv, Misses_oripv + Hits_oripv);
[HR_oripa, ci95_HR_oripa] = binofit(Hits_oripa, Misses_oripa + Hits_oripa);

[HR_oripvv, ci95_HR_oripvv] = binofit(Hits_oripvv, Misses_oripvv + Hits_oripvv);
[HR_oripaa, ci95_HR_oripaa] = binofit(Hits_oripaa, Misses_oripaa + Hits_oripaa);

[FR_oricpv, ci95_FR_oricpv] = binofit(Hits_oricpv, Misses_oricpv + Hits_oricpv);
[FR_oricpa, ci95_FR_oricpa] = binofit(Hits_oricpa, Misses_oricpa + Hits_oricpa);

[FR_oricpvv, ci95_FR_oricpvv] = binofit(Hits_oricpvv, Misses_oricpvv + Hits_oricpvv);
[FR_oricpaa, ci95_FR_oricpaa] = binofit(Hits_oricpaa, Misses_oricpaa + Hits_oricpaa);

n_ori = [Hits_v+Misses_v; Hit_av+Misses_av];
    
end