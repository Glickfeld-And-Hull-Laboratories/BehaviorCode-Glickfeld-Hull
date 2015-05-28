addpath('C:\Users\andrew\Documents\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis');
% Data Path Formatting
laserPowerCell = {};
outcomeCell = {};

disp('Loading data from MWorks Data folder...')
dataPath = 'C:\Users\andrew\Desktop\LaserAnalysis';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
    
% Next line skips ".", "..", and ".DS_Store" hidden files. May not be applicable
% in systems where hidden files are still hidden and should be removed if
% it leads to data skippage.
desNames = fileNames(3:end);
nNames = length(desNames);
for iN = 1:nNames;
        tName = fullfile(dataPath, desNames{iN});
        disp(tName)
        ds=mwLoadData(tName, 'max');
        dailyOutcomes = ds.trialOutcomeCell;
        outcomeCell = [outcomeCell dailyOutcomes];
        laserPower = ds.tTrialLaserPowerMw;
        laserPowerCell = [laserPowerCell laserPower];
        disp('Daily Data Loaded')  
end
%%
disp('Loading complete!')
[sortedPowers, ix] = sort(cell2mat_padded(laserPowerCell));
sortedOutcomes = outcomeCell(ix);
successIx = strcmp(sortedOutcomes, 'success');
earlyIx = strcmp(sortedOutcomes, 'failure');
lateIx = strcmp(sortedOutcomes, 'ignore');
powerLevels = unique(sortedPowers);
powerLevels = powerLevels(1:end-1)
nPowers = length(powerLevels);
percentCorrect = zeros(1,nPowers);
percentEarly = zeros(1,nPowers);
percentLate = zeros(1,nPowers);
for j = 1:nPowers,
    powerIx = sortedPowers == powerLevels(j);
    nTrials = sum(powerIx);
    percentCorrect(j) = sum((powerIx & successIx))/nTrials;
    percentEarly(j) = sum((powerIx & earlyIx))/nTrials;
    percentLate(j) = sum((powerIx & lateIx))/nTrials;
end
hold on
plot(powerLevels, percentLate, 'm-v', 'LineWidth', 2)
plot(powerLevels, percentCorrect, 'k-v', 'LineWidth', 2)
plot(powerLevels, percentEarly, 'c-v', 'LineWidth', 2)
grid on
title('Effect of Laser Stimulation Intensity on Behavioral Trial Outcome')
xlabel('Laser Power (mW)')
ylabel('Probability of Outcome')
clipclip