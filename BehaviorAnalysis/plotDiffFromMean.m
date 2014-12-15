ds = input;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
[num, txt, raw] = xlsread('~/Desktop/LidocaineBook.xlsx');
dates = num(:,1);
date = str2num(ds.saveTime(1:6));
rw = find(dates==date);
trialVals = raw(rw+1,:);
hts = double(cell2mat_padded(input.holdTimesMs));

if ~isnan(trialVals{2}) & ~isnan(trialVals{3}),
    cIx = trialVals{2}:trialVals{3};
    mu = nanmean(hts(cIx));
    cntHolds = hts(cIx)/mu;
    scatter(cntHolds, 'k')
    hold on
end

if ~isnan(trialVals{4}) & ~isnan(trialVals{5}),
    s1Ix = trialVals{4}:trialVals{5};
    s1Holds = hts(s1Ix)/mu;
    scatter(s1Holds, 'b')
    hold on
end

if ~isnan(trialVals{6}) & ~isnan(trialVals{7}),
    s2Ix = trialVals(6):trialVals(7);
    s2Holds = hts(s2Ix)/mu;
    scatter(s2Holds, 'Color', [0 0 0.5])
    hold on
end

if ~isnan(trialVals{8}) & ~isnan(trialVals{9}),
    l1Ix = trialVals{8}:trialVals{9};
    l1Holds = hts(l1Ix)/mu;
    scatter(l1Holds, 'g')
    hold on
end

if ~isnan(trialVals{10}) & ~isnan(trialVals{11}),
    l2Ix = trialVals{10}:trialVals{11};
    l2Holds = hts(l2Ix)/mu;
    plot(l2Holds, 'Color', [0 0.5 0])
    hold on
end