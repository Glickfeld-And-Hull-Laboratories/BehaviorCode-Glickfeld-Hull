addpath('C:\Users\andrew\Documents\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis');
% Data Path Formatting
global dataStru
dataStru = struct;
holds = {};
trs = {};
tooFasts = {};

disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\andrew\Desktop\fixedData\500';
%% 
dirListing = dir(dataPath);
fileNames= {dirListing.name};
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
fileNames = fileNames(3:end);

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
desNames = fileNames;
nNames = length(desNames);
for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        disp(tName) 
        ds=mwLoadData(tName, 'max');
        
        hts = double(cell2mat(ds.holdTimesMs));
        shortIx = hts<200;
        %comment to ignore holds less than 200ms
        hts(hts<200) = NaN;
        
        req = double(cell2mat(ds.tTotalReqHoldTimeMs));
        %normalizedHolds = hts./req
        normalizedHolds = (ds.reactTimesMs);
        tooFast = (ds.tooFastTimeMs);%+req)./req;        
        
        trialNumbers = 1:length(hts);
        
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        normalizedHolds = normalizedHolds(corrIx);
        trialNumbers = trialNumbers(corrIx);
        %tooFast = tooFast(corrIx);
    
        holds = [holds normalizedHolds];
        trs = [trs trialNumbers];
        %tooFasts = [tooFasts tooFast];
        tooFasts = double(ds.tooFastTimeMs);
        
       disp('Daily Data Loaded')  
end

beep
disp('Loading complete!')

holds = cell2mat(holds);
trs = cell2mat(trs);
%tooFasts = cell2mat(tooFasts);
[sortedTrs, ix] = sort(trs);
sortedHolds = holds(ix);
%sortedTooFasts = tooFasts(ix);

%%
holdMeans = {};
holdSTDs = {};
x = 0:10:1990;
y = 10:10:2000;
i=1;
for i = 1:length(x),
    
    iX = x(i)<sortedTrs & sortedTrs<y(i);
    trIx = sortedHolds(iX);
    holdSTDs{i} = nanstd(double(trIx));
    holdMeans{i} = nanmean(double(trIx));
end
hold on
%plot(sortedTooFasts, 'r')
plot([0 1200], [tooFasts tooFasts], 'r');
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
grid on
%plot([0 1200],[1 1], 'k')
xlim([0 1200])
%ylabel('Normalized Hold Time (Actual Hold Time / Required Hold Time)')
ylabel('Mean Reaction Time (ms)')
ylim([0 750])
title('Fixed Time Summary (Corrects Only) -- 500 ms')
xlabel('Trials')