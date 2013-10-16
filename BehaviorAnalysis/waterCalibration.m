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
        fName = strcat('~/Documents/MWorks/Data/data-', char(waterStru(x).values.subjN(i)), '-', mat2str(waterStru(x).values.date(i)), '.mat')
        ds = mwLoadData(fName, 'max');
        waterStru(x).values.avgJuiceTime(i,1) = mean(cell2mat(ds.juiceTimesMsCell));
    end
    waterStru(x).values.waterPerCorrect = waterStru(x).values.water./waterStru(x).values.nCorrect;
end
beep

%%

R1Cell = cat(2, waterStru(1).values.avgJuiceTime, waterStru(1).values.waterPerCorrect);
R1Cell = sortrows(R1Cell);
R2Cell = cat(2, waterStru(2).values.avgJuiceTime, waterStru(2).values.waterPerCorrect);
R2Cell = sortrows(R2Cell);
R3Cell = cat(2, waterStru(3).values.avgJuiceTime, waterStru(3).values.waterPerCorrect);
R3Cell = sortrows(R3Cell);
R4Cell = cat(2, waterStru(4).values.avgJuiceTime, waterStru(4).values.waterPerCorrect);
R4Cell = sortrows(R4Cell);

hold on
r1p = plot(R1Cell(:,1), R1Cell(:,2), 'r')
r2p = plot(R2Cell(:,1), R2Cell(:,2), 'g')
r3p = plot(R3Cell(:,1), R3Cell(:,2), 'b')
r4p = plot(R4Cell(:,1), R4Cell(:,2), 'k')
axis tight

