function [x] = plotWeightWaterHADC8(subjs);
%This function graphs weight and water for the HADC8 Mice
%Written by Andrew McKinney 140131
[numSheet, strSheet, allSheet] = xlsread('~/Downloads/GlickfeldLabTrainingSpreadsheet.xlsx');
global dataStru
%%Quick XLSX sheet casting
xlsSubj = allSheet(:,1);
xlsSubj = xlsSubj(2:end, :);
xlsDate = numSheet(:,1);
xlsEarned = numSheet(:,9);
xlsTotal = numSheet(:,11);
xlsWeight = numSheet(:,2);

for i=1:length(subjMat),   
    %% Switch over subject number to iXXX format if necessary.
    a = subjMat(i);
    if a == [502]
        subStr = 'i502';
    elseif a == [503]
        subStr = 'i503';
    elseif a == [504]
        subStr = 'i504';
    elseif a == [505]
        subStr = 'i505';
    end
    %% Clear Weight/Earned/Total in dataStru for incoming data.
    dataStru(a).values.weight = {};
    dataStru(a).values.earnedWaterMl = {};
    dataStru(a).values.totalWaterMl = {};
    
    %% Subject indexing
    
    outR = regexp(dataStru(a).check.downloadedname, 'data-i([0-9]*)-([0-9]*).mat', 'tokens');
    outRV = cat(1, outR{:});
    outRM = cat(1, outRV{:});
    dateStrList = str2double(outRM(:,2));
    
    subjIx = strcmp(xlsSubj, subStr);
    nDates = length(dateStrList);
    
    subjDates = xlsDate(subjIx);
    subjEarned = xlsEarned(subjIx);
    subjTotal = xlsTotal(subjIx);
    subjWeight = xlsWeight(subjIx);
    
    %% Find weight and earned/total water values based on dataStru dates.
    for ii=1:nDates,
        dataRow = find(subjDates==dateStrList(ii));
        if isempty(dataRow)
            dataStru(a).values.weight{ii} = [NaN];
            dataStru(a).values.earnedWaterMl{ii} = [NaN];
            dataStru(a).values.totalWaterMl{ii} = [NaN];
        elseif ~isempty(dataRow)
            dataStru(a).values.weight{ii} = subjWeight(dataRow);
            dataStru(a).values.earnedWaterMl{ii} = subjEarned(dataRow);
            dataStru(a).values.totalWaterMl{ii} = subjTotal(dataRow);
        end
    end
    dataStru(a).values.weight = cell2mat(dataStru(a).values.weight);
    dataStru(a).values.earnedWaterMl = cell2mat(dataStru(a).values.earnedWaterMl);
    dataStru(a).values.totalWaterMl = cell2mat(dataStru(a).values.totalWaterMl);
end