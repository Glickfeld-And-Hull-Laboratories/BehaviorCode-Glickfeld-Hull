function [regPlot targRewMs] = rewardCalculator(nRig, nCorrect, earnedRewMl, minCorrect)
%% rewardCalculator
%   Given [regPlot targRewMs] = rewardCalculator(nRig, nCorr, earnedRewMl)
% nRig = Rig number (1-4)
% nCorr = Expected number of correct trials
% earnedRewMl = Expected water for mouse to earn in milliliters
% minCorrect = Minimum daily number of corrects. Days with fewer than 
%minCorrect corrects will be ignored and will not be used in regression.
%
% Function searches much like waterCalibration, finding values for a single
%rig, then using those values to calculate

if nargin < 3,
    error('Insufficient Input Arguments')
elseif nargin == 3,
    % If minCorrect isn't set, will sort out days with less than 50 nCorr.
    minCorrect = 50;
end

%% Indexing/Constants Evaluation
[numSheet, strSheet, allSheet] = xlsread('~/Downloads/HullTrainingSpreadsheet1.xlsx');
rewReg = struct;
xlsSubj = allSheet(2:end,1);
xlsDate = numSheet(:,1);
xlsRig = numSheet(:,4);
xlsWater = numSheet(:,9);
xlsCorrect = numSheet(:, 8);
rigIx = xlsRig == nRig;
rewReg(1).values.water = xlsWater(rigIx);
rewReg(1).values.date = xlsDate(rigIx);
rewReg(1).values.subjN = xlsSubj(rigIx);
rewReg(1).values.nCorrect = xlsCorrect(rigIx);
%% Data Loading and Averaging
for i = 1:length(rewReg(1).values.date),
    fName = strcat('~/Documents/MWorks/Data/data-', char(rewReg(1).values.subjN(i)), '-', mat2str(rewReg(1).values.date(i)), '.mat');
    ds = mwLoadData(fName, 'max');
    rewReg(1).values.avgJuiceTime(i,1) = nanmean(cell2mat(ds.juiceTimesMsCell));
end
rewReg(1).values.waterPerCorrect = rewReg(1).values.water./rewReg(1).values.nCorrect;
%% Sorting into 10ms bins, averaging, and throwing out nCorr < minCorrect
RCell = cat(2, rewReg(1).values.avgJuiceTime, rewReg(1).values.waterPerCorrect, rewReg(1).values.nCorrect);
RCell = sortrows(RCell);
ltIx = RCell(:,3)<minCorrect;
RJuice = RCell(:,1);
RJuice = RJuice(~ltIx);
RJuiceA = roundn(RJuice, 1);
uIx = unique(RJuiceA);
R1WperC = RCell(:,2);
R1WperC = R1WperC(~ltIx);

for i=1:length(uIx),
    AvgJ(i) = nanmean(R1WperC(find(RJuiceA==uIx((i)))));
end

%% Plotting
rigPParamCell = {'r.-' 'g.-' 'b.-' 'k.-'};
rigPParam = rigPParamCell{nRig};
hold on
rigP = plot(uIx, AvgJ.*1000, rigPParam);
todayStr = datestr(today, 'dd mmmm yyyy');
tName = strcat('Reward Regression Plot - Rig', num2str(nRig), ' - Generated:  ', todayStr);
title(tName);
xlabel('Single Reward Size (ms)')
ylabel('Saccharin Solution Dispensed (\muL)')
%% 3rd Order Polynomial Regression
avgJ = (AvgJ.*1000)';
p = polyfit(uIx, avgJ, 3);
f = polyval(p, uIx);
regPParamCell = {'r-.' 'g-.' 'b-.' 'k-.'};
regPParam = regPParamCell{nRig};
hold on
regPlot = plot(uIx, f, regPParam);
axis tight
grid on
legend('Calibration Curve', 'Regression Curve')
sName = strcat('rewardCalculator-R', num2str(nRig), '-', datestr(today, 'yymmdd'), '.pdf');
    saveas(gcf, sName);
    beep
%% Evaluation of polynomial regression to yield targRewMs value

eMeanSingleRewUl = (earnedRewMl/nCorrect).*1000;
p2 = polyfit(avgJ, uIx, 3);
targRewMs = polyval(p2, eMeanSingleRewUl);





