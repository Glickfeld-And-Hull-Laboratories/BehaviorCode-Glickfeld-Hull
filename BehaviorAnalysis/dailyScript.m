function [ output_args ] = dailyScript(subjMat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
startup

global dataStru
plotWeightWaterHADC8([ 502 503 506 507 510 511]);
dataStru = struct;
dataStru = dataLoad(subjMat);
dataStru = wwLoad(subjMat);
plotPerformance(subjMat);
quit
end

