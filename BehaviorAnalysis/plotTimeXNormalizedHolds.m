holds = {};
starts = {};
for i=1:length(ii),
    x = ii{i};
    
    hts = double(cell2mat(x.holdTimesMs));
    req = double(cell2mat(x.tTotalReqHoldTimeMs));
    normalizedHolds = hts./req;
    
    htStarts = x.holdStartsMs;
    holdStarts = (cell2mat(htStarts) - htStarts{1})/60000;
    
    holds = [holds normalizedHolds];
    starts = [starts holdStarts];
end
holds = cell2mat(holds);
starts = cell2mat(starts);

[sortedStarts, ix] = sort(starts);
sortedHolds = holds(ix);

holdMeans = {};
holdSTDs = {};
x = 0:2:100;
y = 2:2:102;
i=1;
for i = 1:length(x),
    
    iX = x(i)<sortedStarts & sortedStarts<y(i)
    trs = sortedHolds(iX);
    holdSTDs{i} = std(trs);
    holdMeans{i} = nanmean(trs);
end
hold on
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'g-x')
grid on
ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)')
plot([0 100],[1 1])
plot([0 100],[1 1], 'k')
ylim([0.75 1.25])
title('i226 - Binned Data from 150302 to 150305')
clipclip
xlabel('Minutes in Training')