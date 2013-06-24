function [ output ] = exptPlotPdfOffline(subjNum, dateString, dataBlock, ...
                                         plotFunctionFullPath, doSavePdf)
% exptPlotPdfOffline
% 
%       Draw pdf figure from that day
%       plotFunctionFullPath: path to plotting function, default to ~/ExperimentXML-git/HoldAndDetectConstant/MATLAB_window_files/plotOnlineHist.m
%
% 

% histed 130116, original by amckinney

if nargin < 3, dataBlock = 'max'; end
if nargin < 4, plotFunctionFullPath = []; end
if nargin < 5, doSavePdf = false; end

cs = exptConstants;
addpath('~/Repositories/MWorksMatlabToolbox');
addpath('~/Repositories/tools-mh');

if isempty(plotFunctionFullPath)
    plotFunctionFullPath = fullfile(cs.homeDir, 'ExperimentXML-git/HoldAndDetectConstant/MATLAB_window_files/plotOnlineHist.m');
end

[tPath tName tExt] = fileparts(plotFunctionFullPath);
addpath(tPath, '-begin');
plotFH = str2func(tName);

%dataName= strcat('data-', subjNum, '-', dateStrng, '.mat');
dataName = sprintf('data-i%03d-%s.mat', subjNum, dateString);
fileName= fullfile(cs.dataPath, dataName);

ds=mwLoadData(fileName, dataBlock, true);

figH=figure(4);
set(figH, 'Name', mfilename);

feval(plotFH, ds, ds);

bH = findobj(gcf, 'Type', 'uicontrol');
assert(~isempty(bH));
delete(bH)

if doSavePdf
    outName = '~/Desktop/plotPdfOffline.pdf';
    fprintf(1,'Saving PDF file to %s...', outName);

    exportfig_print(gcf, outName, ...
                 'FileFormat', 'pdf', ...
                 'Size', [12 12], ...
                 'PrintUI', false);
    fprintf(1, ' done.\n');
end

rmpath(tPath);