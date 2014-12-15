%function [pdf] = plotLidocaineTDiff(input),
% Given an input data structure, will
ds = input;
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
[num, txt, raw] = xlsread('~/Desktop/LidocaineBook.xlsx');
dates = num(:,1);
date = str2num(ds.saveTime(1:6));
rw = find(dates==date);
trialVals = raw(rw+1,:);


thisT = double(cell2mat_padded(ds.tThisTrialStartTimeMs));
%thisT = thisT(2:end);
lastT = double(cell2mat_padded(ds.tLastTrialStartTimeMs));
%lastT = lastT(2:end);
diffT = thisT-lastT;

if ~isnan(trialVals{2}) & ~isnan(trialVals{3}),
    cIx = trialVals{2}:trialVals{3};
    cntHolds = diffT(cIx);
    plot(cntHolds(2:end), 'k')
    hold on
end

if ~isnan(trialVals{4}) & ~isnan(trialVals{5}),
    s1Ix = trialVals{4}:trialVals{5};
    s1Holds = diffT(s1Ix);
    plot(s1Holds(2:end), 'b')
    hold on
end

if ~isnan(trialVals{6}) & ~isnan(trialVals{7}),
    s2Ix = trialVals{6}:trialVals{7};
    s2Holds = diffT(s2Ix);
    plot(s2Holds(2:end), 'g')
    hold on
end

if ~isnan(trialVals{8}) & ~isnan(trialVals{9}),
    l1Ix = trialVals{8}:trialVals{9};
    l1Holds = diffT(l1Ix);
    plot(l1Holds(2:end), 'COlor', [1 0.6 0])
    hold on
end

if ~isnan(trialVals{10}) & ~isnan(trialVals{11}),
    l2Ix = trialVals{10}:trialVals{11};
    l2Holds = diffT(l2Ix);
    plot(l2Holds(2:end), 'r')
    hold on
end
grid on
xlabel('Trial Number')
ylabel('Intertrial Start Time Difference (ms)')
outPath = '~/Desktop/';
title(ds.savedDataName(40:end))
fName = strcat('trialStartDiff-', ds.savedDataName(40:55), '.pdf')
fileName = strcat(outPath, fName);
exportfig_print(gcf, fileName, 'FileFormat', 'pdf', ...
             'Size', [12 9], ...
             'PrintUI', false)
pdf = gcf;
