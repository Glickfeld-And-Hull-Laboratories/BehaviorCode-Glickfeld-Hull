if ispc,
    disp('I AM A PC')
    addpath('C:\Users\jake\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis\');
elseif isunix,
    disp('I AM A MAC OR LINUX COMPUTER')
    addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
end
% Data Path Formatting
global dataStru
dataStru = struct;
holds = {};
starts = {};
tooFasts = {};

disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\jake\Desktop\AnalysisFolder\temp';
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
        hts(hts<200) = NaN;
        req = double(cell2mat(ds.tTotalReqHoldTimeMs));
        normalizedHolds = (ds.reactTimesMs);%hts./req;
        tooFasts = ((ds.tooFastTimeMs));%+req)./req;
        
    
        htStarts = ds.holdStartsMs;
        holdStarts = (cell2mat(htStarts) - htStarts{1})/60000;
        
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        normalizedHolds = normalizedHolds(corrIx);
        holdStarts = holdStarts(corrIx);
        %tooFast = tooFast(corrIx);
    
        holds = [holds normalizedHolds];
        starts = [starts holdStarts];
        %tooFasts = [tooFasts tooFast];
        
       disp('Daily Data Loaded')  
end
%%
beep
disp('Loading complete!')

holds = cell2mat(holds);
starts = cell2mat(starts);
%tooFasts = cell2mat(tooFasts);

[sortedStarts, ix] = sort(starts);
sortedHolds = holds(ix);
%sortedTooFasts = tooFasts(ix);

%%
holdMeans = {};
holdSTDs = {};
x = 0:2:118;
y = 2:2:120;
for i = 1:length(x),
    
    iX = x(i)<sortedStarts & sortedStarts<y(i);
    trs = sortedHolds(iX);
    holdSTDs{i} = nanstd(double(trs));
    holdMeans{i} = nanmean(double(trs));
end
hold on
%plot(sortedTooFasts, 'r')
%plot([0 120],[1 1], 'k')
plot([0 1200], [tooFasts tooFasts], 'r');
errorbar(y, cell2mat(holdMeans), cell2mat(holdSTDs), 'k-x')
grid on
ylabel('Mean Hold Time (Actual Hold Time / Required Hold Time)');
%ylabel('Mean Reaction Time (ms)')
ylim([0 500])
xlim([1 100])
title('Fixed Time Summary (Corrects Only) -- 500 ms ')
xlabel('Minutes in Training')