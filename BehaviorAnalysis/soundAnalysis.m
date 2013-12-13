function [soundStru] = soundAnalysis(subj, date)
%% Function to analyze performance at different delaySoundVolumes in a day

if nargin<2,
    date = datestr(today, 'yymmdd');
end

dataPath = '~/Documents/MWorks/Data';
fName = strcat(dataPath, '/data-i', subj, '-', date, '.mat');
ds = mwLoadData(fName, 'max');


vol(1) = ds.firstTrConsts.delaySoundVolume
ix = regexp(ds.changedStr, 'Trial *([0-9]*)', 'tokens')

end
