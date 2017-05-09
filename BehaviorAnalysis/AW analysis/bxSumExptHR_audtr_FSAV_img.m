if hr_calc == 1
[h_amp, edge, bin_amp] = histcounts(audTargets, amp_edges);
% [h_oric, edge, bin_oric] = histcounts(invVisTargets, amp_edges);

amp_bins = unique(bin_amp);
avg_amp = zeros(1,max(amp_bins,[],2));
sem_amp = zeros(1,max(amp_bins,[],2));
Hits_amp = zeros(1,max(amp_bins,[],2));
Misses_amp = zeros(1,max(amp_bins,[],2));
% avg_oric = zeros(1,max(amp_bins,[],2));
% sem_oric = zeros(1,max(amp_bins,[],2));
% FAs_ori = zeros(1,max(amp_bins,[],2));
% CRs_ori = zeros(1,max(amp_bins,[],2));
rct_amp = zeros(1,max(amp_bins,[],2));
rctsem_amp = zeros(1,max(amp_bins,[],2));
% rct_oric = zeros(1,max(amp_bins,[],2));
% rctsem_oric = zeros(1,max(amp_bins,[],2));
n_amp = zeros(2,max(amp_bins,[],2));
for ibin = 1:max(amp_bins,[],2)
    ind = find(bin_amp == ibin);
    if length(ind)>=1
        avg_amp(ibin) = mean(audTargets(:,ind),2);
        sem_amp(ibin) = std(audTargets(:,ind),[],2)./sqrt(length(ind));
        Hits_amp(ibin) = sum(successIx(ind),2);
        Misses_amp(ibin) = sum(missedIx(ind),2);
%             rct_amp(ibin) = mean(targetReact(successIx(ind)),2)
%             rctsem_amp(ibin) = std(targetReact(successIx(ind)),[],2)./length(ind);
    end
% %     indc = find(bin_oric == ibin);
% %     if length(indc)>=1
% %         avg_oric(ibin) = mean(invVisTargets(:,indc),2);
% %         sem_oric(ibin) = std(invVisTargets(:,indc),[],2)./sqrt(length(indc));
% %         FAs_ori(ibin) = sum(invHitIx(indc),2);
% %         CRs_ori(ibin) = sum(invMissIx(indc),2);
% % %             rct_oric(ibin) = mean(catchReact(invHitIx(indc)),2);
% % %             rctsem_oric(ibin) = std(catchReact(invHitIx(indc)),[],2)./length(indc);
% %     end
end
[HR_amp, ci95_HR_amp] = binofit(Hits_amp, Misses_amp + Hits_amp);
% [FR_ori, ci95_FR_ori] = binofit(FAs_ori, CRs_ori + FAs_ori);
n_amp = [Hits_amp+Misses_amp];

%%
elseif hr_calc == 2
    
  
[h_amp, edge, bin_amp] = histcounts(audTargets_all, amp_edges);
% [h_oric, edge, bin_oric] = histcounts(invaudTargets_all, amp_edges);

amp_bins = unique(bin_amp);
avg_amp = zeros(1,max(amp_bins,[],2));
sem_amp = zeros(1,max(amp_bins,[],2));
Hits_amp = zeros(1,max(amp_bins,[],2));
Misses_amp = zeros(1,max(amp_bins,[],2));
% avg_oric = zeros(1,max(amp_bins,[],2));
% sem_oric = zeros(1,max(amp_bins,[],2));
% FAs_ori = zeros(1,max(amp_bins,[],2));
% CRs_ori = zeros(1,max(amp_bins,[],2));
rct_amp = zeros(1,max(amp_bins,[],2));
rctsem_amp = zeros(1,max(amp_bins,[],2));
% rct_oric = zeros(1,max(amp_bins,[],2));
% rctsem_oric = zeros(1,max(amp_bins,[],2));
n_amp = zeros(2,max(amp_bins,[],2));
for ibin = 1:max(amp_bins,[],2)
    ind = find(bin_amp == ibin);
    if length(ind)>=1
        avg_amp(ibin) = mean(audTargets_all(:,ind),2);
        sem_amp(ibin) = std(audTargets_all(:,ind),[],2)./sqrt(length(ind));
        Hits_amp(ibin) = sum(successIx_all(ind),2);
        Misses_amp(ibin) = sum(missedIx_all(ind),2);
%             rct_amp(ibin) = mean(targetReact(successIx(ind)),2)
%             rctsem_amp(ibin) = std(targetReact(successIx(ind)),[],2)./length(ind);
    end
% %     indc = find(bin_oric == ibin);
% %     if length(indc)>=1
% %         avg_oric(ibin) = mean(invaudTargets_all(:,indc),2);
% %         sem_oric(ibin) = std(invaudTargets_all(:,indc),[],2)./sqrt(length(indc));
% %         FAs_ori(ibin) = sum(invHitIx_all(indc),2);
% %         CRs_ori(ibin) = sum(invMissIx_all(indc),2);
% % %             rct_oric(ibin) = mean(catchReact(invHitIx(indc)),2);
% % %             rctsem_oric(ibin) = std(catchReact(invHitIx(indc)),[],2)./length(indc);
% %     end
end
[HR_amp, ci95_HR_amp] = binofit(Hits_amp, Misses_amp + Hits_amp);
% [FR_ori, ci95_FR_ori] = binofit(FAs_ori, CRs_ori + FAs_ori);
n_amp = [Hits_amp+Misses_amp];
    
end