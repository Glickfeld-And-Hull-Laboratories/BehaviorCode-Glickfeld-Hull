[numSheet, strSheet, allSheet] = xlsread('~/Downloads/HullMouseTrainingProtocol.xlsx');
xlsSubj = allSheet(:,1);
xlsSubj = xlsSubj(2:end, :);
xlsDate = numSheet(:,1);
xlsRig = numSheet(:,4);
xlsWater = numSheet(:,9);
waterStru = struct;

r1Ix = xlsRig == 1;
r2Ix = xlsRig == 2;
r3Ix = xlsRig == 3;
r4Ix = xlsRig == 4;

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


for x = 1:4,
    for i = 1:length(waterStru(x).values.date,
        fName = strcat('~/Documents/MWorks/Data/data-', waterStru(x).values.subjN(i), '-', waterStru(x).values.date(i), '.mat')
        ds = mwLoadData(fName, 'max');
        waterStru(x).values.avgJuiceTime = mean(cell2mat(ds.juiceTimesMsCell));
    end
end
