function [ pdf ] = fitNormalToHolds(ds)
%FITNORMALTOHOLDS: reads the hold times of a single day's RewardDelay data
%and fits a normal distribution and plots relative to the onset of the
%reward window.
if isfield(ds, 'tDelayTimeMs'),
    expID = 'RewardDelay';
end
holdTimes = cell2mat_padded(ds.holdTimesMs);
ht = sort(double(holdTimes));
htMax = ceil(max(ht)/100)*100;
bins = 0:100:htMax;

[a] = histfit(ht, 50);
mu = mean(a(2).XData)
end
