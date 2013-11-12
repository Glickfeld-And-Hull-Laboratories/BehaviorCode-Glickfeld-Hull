function [ dataStru ] = wwLoad(subjMat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%subjMat = [1 5 6 7 202 205 206 209];
[numSheet, strSheet, allSheet] = xlsread('~/Downloads/HullMouseTrainingProtocol.xlsx');
global wwStru
wwStru = struct;
global dataStru
%%Quick XLSX sheet casting
xlsSubj = allSheet(:,1);
xlsSubj = xlsSubj(2:end, :);
xlsDate = numSheet(:,1);
xlsEarned = numSheet(:,9);
xlsTotal = numSheet(:,11);
xlsWeight = numSheet(:,2);

for i=1:length(subjMat),   
    %% Quick bit to switch over subject number to iXXX format.
    a = subjMat(i);
    if a == [1]
        subStr = 'i001';
    elseif a == [2]
        subStr = 'i002';
    elseif a == [3]
        subStr = 'i003';
    elseif a == [5]
        subStr = 'i005';
    elseif a == [6]
        subStr = 'i006';
    elseif a == [7]
        subStr = 'i007';
    elseif a == [202]
        subStr = 'i202';
    elseif a == [203]
        subStr = 'i203';
    elseif a == [205]
        subStr = 'i205';
    elseif a == [206]
        subStr = 'i206';
    elseif a == [209]
        subStr = 'i209';
    end
    
    %% Subject indexing
    subjIx = strcmp(xlsSubj, subStr);
    nDates = length(dataStru(a).dates);
    subjDates = xlsDate(subjIx);
    subjEarned = xlsEarned(subjIx);
    subjTotal = xlsTotal(subjIx);
    subjWeight = xlsWeight(subjIx);
    
    %% Find weight and earned/total water values based on dataStru dates.
    for ii=1:nDates,
        dataRow = find(subjDates==dataStru(a).dates(ii));
        if isempty(dataRow)
            dataStru(a).values.weight{ii} = [];
            dataStru(a).values.earnedWaterMl{ii} = [];
            dataStru(a).values.totalWaterMl{ii} = [];
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

