%%
[numSheet, strSheet, allSheet] = xlsread('~/Downloads/HullMouseTrainingProtocol.xlsx');
waterStru = struct;
%%
xlsSubj = allSheet(:,1);
xlsSubj = xlsSubj(2:end, :);
xlsDate = numSheet(:,1);
xlsRig = numSheet(:,4);
xlsWater = numSheet(:,9);
xlsCorrect = numSheet(:, 8);

%%
r1Ix = xlsRig == 1;
r2Ix = xlsRig == 2;
r3Ix = xlsRig == 3;
r4Ix = xlsRig == 4;
%%
waterStru(1).values.water = xlsWater(r1Ix);
waterStru(2).values.water = xlsWater(r2Ix);
waterStru(3).values.water = xlsWater(r3Ix);
waterStru(4).values.water = xlsWater(r4Ix);

waterStru(1).values.date = xlsDate(r1Ix);
waterStru(2).values.date = xlsDate(r2Ix);
waterStru(3).values.date = xlsDate(r3Ix);
waterStru(4).values.date = xlsDate(r4Ix);

waterStru(1).values.subjN = xlsSubj(r1Ix);
waterStru(2).values.subjN = xlsSubj(r2Ix);
waterStru(3).values.subjN = xlsSubj(r3Ix);
waterStru(4).values.subjN = xlsSubj(r4Ix);

waterStru(1).values.nCorrect = xlsCorrect(r1Ix);
waterStru(2).values.nCorrect = xlsCorrect(r2Ix);
waterStru(3).values.nCorrect = xlsCorrect(r3Ix);
waterStru(4).values.nCorrect = xlsCorrect(r4Ix);

%%
for x = 1:4,
    for i = 1:length(waterStru(x).values.date),
        fName = strcat('~/Documents/MWorks/Data/data-', char(waterStru(x).values.subjN(i)), '-', mat2str(waterStru(x).values.date(i)), '.mat');
        ds = mwLoadData(fName, 'max');
        waterStru(x).values.avgJuiceTime(i,1) = mean(cell2mat(ds.juiceTimesMsCell));
    end
    waterStru(x).values.waterPerCorrect = waterStru(x).values.water./waterStru(x).values.nCorrect;
end
beep

%%

R1Cell = cat(2, waterStru(1).values.avgJuiceTime, waterStru(1).values.waterPerCorrect, waterStru(1).values.nCorrect);
R1Cell = sortrows(R1Cell);
ltIx1 = R1Cell(:,3)<75;
R1Juice = R1Cell(:,1);
R1Juice = R1Juice(~ltIx1);

R1JuiceA = roundn(R1Juice, 1);
uIx1 = unique(R1JuiceA);

R1WperC = R1Cell(:,2);
R1WperC = R1WperC(~ltIx1);

for i=1:length(uIx1),
    AvgJ1(i) = nanmean(R1WperC(find(R1JuiceA==uIx1((i)))));
end

R2Cell = cat(2, waterStru(2).values.avgJuiceTime, waterStru(2).values.waterPerCorrect, waterStru(2).values.nCorrect);
R2Cell = sortrows(R2Cell);
ltIx2 = R2Cell(:,3)<75;
R2Juice = R2Cell(:,1);
R2Juice = R2Juice(~ltIx2);
R2JuiceA = roundn(R2Juice, 1);
uIx2 = unique(R2JuiceA);
R2WperC = R2Cell(:,2);
R2WperC = R2WperC(~ltIx2);

for i=1:length(uIx2),
    AvgJ2(i) = nanmean(R2WperC(find(R2JuiceA==uIx2((i)))));
end


R3Cell = cat(2, waterStru(3).values.avgJuiceTime, waterStru(3).values.waterPerCorrect, waterStru(3).values.nCorrect);
R3Cell = sortrows(R3Cell);
ltIx3 = R3Cell(:,3)<75;
R3Juice = R3Cell(:,1);
R3Juice = R3Juice(~ltIx3);
R3JuiceA = roundn(R3Juice, 1);
uIx3 = unique(R3JuiceA);
R3WperC = R3Cell(:,2);
R3WperC = R3WperC(~ltIx3);

for i=1:length(uIx3),
    AvgJ3(i) = nanmean(R3WperC(find(R3JuiceA==uIx3((i)))));
end

R4Cell = cat(2, waterStru(4).values.avgJuiceTime, waterStru(4).values.waterPerCorrect, waterStru(4).values.nCorrect);
R4Cell = sortrows(R4Cell);
ltIx4 = R4Cell(:,3)<75;
R4Juice = R4Cell(:,1);
R4Juice = R4Juice(~ltIx4);
R4JuiceA = roundn(R4Juice, 1);
uIx4 = unique(R4JuiceA);
R4WperC = R4Cell(:,2);
R4WperC = R4WperC(~ltIx4);

for i=1:length(uIx4),
    AvgJ4(i) = nanmean(R4WperC(find(R4JuiceA==uIx4((i)))));
end

hold on
r1p = plot(uIx1, AvgJ1.*1000, 'r.-')
r2p = plot(uIx2, AvgJ2.*1000, 'g.-')
r3p = plot(uIx3, AvgJ3.*1000, 'b.-')
r4p = plot(uIx4, AvgJ4.*1000, 'k.-')
legend('1', '2', '3', '4', 'Location', 'Best')
todayStr = datestr(today, 'dd mmmm yyyy');
tName = strcat('Water Calibration Plot - Generated  ', todayStr);
title(tName);
xlabel('Average Reward Size (ms)')
ylabel('Saccharin Solution Dispensed (\muL)')
axis tight
grid on

sName = strcat('waterCalibration-', datestr(today, 'yymmdd'), '.pdf');
    saveas(gcf, sName)
    

