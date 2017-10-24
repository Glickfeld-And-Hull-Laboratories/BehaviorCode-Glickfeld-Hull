function [h, p, chistat] = binChiHR(hits1,misses1,hits2,misses2)
% is there a significant difference (by chi-square test) between each bin,
% (bonferonni corrected)
nbin = length(hits1);
total1 = hits1+misses1;
total2 = hits2+misses2;
h = nan(1,nbin);
p = nan(1,nbin);
chistat = nan(1,nbin);
for ib = 1:nbin
    X = [hits1(ib) hits2(ib)];
    N = [total1(ib) total2(ib)];
   [h(ib),p(ib),chistat(ib)] = prop_test(X , N, 0.05/nbin);
end

end