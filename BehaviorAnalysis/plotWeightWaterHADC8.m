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
    
    %% Actual Plotting here
     %Ren
    weight = dataStru(a).values.weight;
    earnedmL = dataStru(a).values.earnedWaterMl;
    totalmL =  dataStru(a).values.totalWaterMl;   
        
    %plot(xDays, totalmL, 'b', 'LineWidth', 2)
    [Ax, H1, H2] = plotyy(xDays, earnedmL, xDays, weight);
    set(H1,'LineWidth', 2);
    set(H1,'Color', 'b');
    set(H2,'LineWidth', 2);
    set(H2,'Color', 'k');
    set(Ax(1), 'XLim', [1 max(xDays)]);
    set(Ax(2), 'XLim', [1 max(xDays)]);
    set(Ax(2), 'YLim', [16 25]);
    set(Ax(1), 'YLim', [0 max(earnedmL)]);
    set(Ax(1),'YColor', [0 0 1])
    set(Ax(2),'YColor', [0 0 0])
    xlabel('Training Day');
    set(get(Ax(1),'YLabel'), 'String', 'Earned Water Volume (mL)');
    set(get(Ax(2),'YLabel'), 'String', 'Weight (g)');
    set(Ax(2), 'YTick', [16:1:25]);
    set(Ax(1), 'YTick', [0:0.2:ceil(max(earnedmL)/0.2)*0.2]);
    grid on
    title('Weight and Water')
    fName = strcat('Performance -- i', num2str(a), '-- Generated:', datestr(today, 'dd mmmm yyyy'));
    set(gcf, 'Name', fName);
    
end