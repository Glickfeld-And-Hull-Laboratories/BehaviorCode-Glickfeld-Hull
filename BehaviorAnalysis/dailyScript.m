function [ output_args ] = dailyScript(subjMat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
startup

global dataStru
dataStru = struct;
dataStru = dataLoad(subjMat);
dataStru = wwLoad(subjMat);
plotPerformance(subjMat);
quit
end

