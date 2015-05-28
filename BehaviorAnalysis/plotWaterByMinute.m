addpath('C:\Users\andrew\Documents\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis');
% Data Path Formatting
global dataStru
dataStru = struct;
totalRewardPercents = {};
totalStarts = {};
totalTotals = {};

dataPath = 'C:\Users\andrew\Desktop\fixedData\2000';
disp('Loading data from MWorks Data folder...')
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
        %format rewards as percent of total reward
        rewards = cell2mat_padded(ds.juiceTimesMsCell);
        rewardPercents = (rewards./sum(rewards));
        totalRewardPercents = [totalRewardPercents rewardPercents];
        totalTotals = [totalTotals sum(rewards)];
        
        trialStarts = ds.holdStartsMs;
        trialStarts = (cell2mat(trialStarts) - trialStarts{1})/60000;
        totalStarts = [totalStarts trialStarts];
        
        subjString = strcat('i', mat2str(ds.subjectNum));
        
        disp('Daily Data Loaded')
end
%%
beep
disp('Loading complete!')

totalRewardPercents = cell2mat_padded(totalRewardPercents)./iN;
%totalRewardPercents(totalRewardPercents == 0) = NaN;
totalStarts = cell2mat_padded(totalStarts);

[sortedStarts, ix] = sort(totalStarts);
sortedRewPercents = totalRewardPercents(ix);

%%
holdMeans = {};
holdSTDs = {};
x = 0:1:119;
y = 1:1:120;
for i = 1:length(x),
    
    iX = x(i)<=sortedStarts & sortedStarts<y(i);
    rews = sortedRewPercents(iX);
    rewSTDs{i} = nanstd(double(rews));
    rewMeans{i} = nanmean(double(rews));
end
%%
hold on
plot(x, cell2mat_padded(rewMeans)*100, 'k-v')
grid on
xlabel('Minutes in Training')
ylabel('Percent of Total Reward Earned')

