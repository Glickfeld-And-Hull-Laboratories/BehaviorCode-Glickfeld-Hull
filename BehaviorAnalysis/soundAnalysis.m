function [vol holdTimes] = soundAnalysis(subj, date)
%% Function to analyze performance at different delaySoundVolumes in a day

if nargin<2,
    date = datestr(today, 'yymmdd');
end

%dataPath = '~/Documents/MWorks/Data';
dataPath = '~/Downloads/i10zip';
fName = strcat(dataPath, '/data-', subj, '-', date, '.mat');
ds = mwLoadData(fName, 'max');

%Set first trial's values to first vol and trChanged values
vol(1) = ds.firstTrConsts.delaySoundVolume;
trChanged(1) = 1;

%Indexes the changedStr cell and pulls out delayVolume change and when it
%occurred
ix3 = regexp(ds.changedStr, 'Trial *([0-9]*): delaySoundVolume: -> (\d*\.?\d*)', 'tokens');
ix2 = cat(1, ix3{:});
ix = cat(1, ix2{:});
nTrs = length(ds.holdTimesMs);

for i=1:length(ix)
    n = i+1;
    thisVol = ix(i,2);
    thisTrN = ix(i,1);
    vol(n,1)= str2double(ix(i,2));
    trChanged(n,1) = str2double(ix(i,1)); 
    holdTimes{i} = cell2mat_padded(ds.holdTimesMs(trChanged(i):trChanged(n)));
    if i==length(ix)
        trChanged(end+1)= nTrs;
        holdTimes{n} = cell2mat_padded(ds.holdTimesMs(trChanged(n):nTrs));
    end
end
for k=1:length(holdTimes)
    gtIx(k) = lt(length(holdTimes{k}), 20);
end
holdTimes = holdTimes(~gtIx);
vol = vol(~gtIx);


% Hold time histogram

if ds.delayTimeMs  < 50,
    maxX = 2000;
    disp('true');
else
    maxFail = double(ds.delayTimeMs+ds.waitForLateTimeMs+(ds.postRewardWindowMs));
    maxX = ceil((maxFail)./500)*500;  % round up to nearest 500 ms.
end

binWidth = 50;	% robust version of std
nBins = ceil(maxX ./ binWidth);

edges = linspace(0, maxX, nBins);

for j=1:length(holdTimes),
    if j==1,
        lColor = 'r';
    elseif j==2,
        lColor = [1 0.64 0]; %orange
    elseif j==3,
        lColor = 'y';
    elseif j==4,
        lColor = 'g';
    elseif j==5,
        lColor = 'b'
    end
    binnedTrs = histc(double(holdTimes{j}), edges);
    bH = plot(edges, binnedTrs, 'Color', lColor, 'LineWidth', 2)
    hold on
end
hold on;
xLim = [0 maxX];
set(gca, 'XLim', xLim);
yLim = get(gca, 'YLim');
if (length(get(gca, 'XTick')) > 4)
  xT = (0:500:maxX);
  set(gca, 'XTick', xT);
end
vH = vert_lines([ds.delayTimeMs-(ds.preRewardWindowMs) ds.delayTimeMs ds.delayTimeMs+(ds.postRewardWindowMs)]);
set(vH, 'Color','k');
axis tight
ylabel('Number of Trials')
xlabel('Lever Hold Time (ms)')
ttl = sprintf('Sound Analysis: Hold Times from %s on %s', subj, date);
title(ttl)

if length(holdTimes)==1,
    legend(num2str(vol(1)));
elseif length(holdTimes)==2,
    legend(num2str(vol(1)), num2str(vol(2)));
elseif length(holdTimes)==3,
    legend(num2str(vol(1)), num2str(vol(2)), num2str(vol(3)));
elseif length(holdTimes)==3,
    legend(num2str(vol(1)), num2str(vol(2)), num2str(vol(3)), num2str(vol(4)));
end

sName = strcat('~/Documents/MATLAB/SoundAnalysis/soundAnalysis-', subj,'-', date, '.pdf');
epParams = { gcf, sName, ...
    'FileFormat', 'pdf', ...
    'Size', [12 12], ...
    'PrintUI', false };
exportfig_print(epParams{:});
%close