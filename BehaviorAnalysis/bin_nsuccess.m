function [hits, misses] = bin_nsuccess(targets,hit_ind,miss_ind,edges,min_trials)

[~,~,bin_ind] = histcounts(targets,edges);

bins = unique(bin_ind);
max_bins = max(bins);
% avg_tar = zeros(1,max_bins);
% sem_tar = zeros(1,max_bins);
hits = zeros(1,max_bins);
misses = zeros(1,max_bins);

for ibin = 1:max_bins
    ind = bin_ind == ibin;
    if length(ind) >= min_trials
        hits(ibin) = sum(hit_ind(ind));
        misses(ibin) = sum(miss_ind(ind));
%         avg_tar(ibin) = mean(targets(ind));
%         sem_tar(ibin) = ste(targets(ind),2);
    end 
end
% nsuccess = hits;
% total = hits+misses;
% [hr, ci95] = binofit(hits,hits+misses);    

end