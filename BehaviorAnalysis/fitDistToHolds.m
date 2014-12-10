function [a b] = fitDistToHolds(ds)
%FITNORMALTOHOLDS: given the hold times of a single day's RewardDelay data, 
% fits a normal distribution and plots relative to the onset of the 
%reward window.

holdTimes = cell2mat_padded(ds.holdTimesMs);
ht = sort(holdTimes);
htMax = ceil(max(ht)/25)*25;
bins = 0:25:ceil(htMax/25)*25;
cnt = histcounts(ht, bins);
[a] = plot(bins(2:end), cnt);
%[a] = histfit(double(ht), length(bins), 'kernel');
%stairs(ht, 'Color', [0 1 0])
%set(a(1), 'FaceColor', [1 0 1])
%set(a(2), 'Visible', 'off')
[b] = vert_lines(ds.delayTimeMs - ds.preRewardWindowMs);
set(b, 'Color', [ 0 1 0])
%mu = mean(a(2).XData)
end