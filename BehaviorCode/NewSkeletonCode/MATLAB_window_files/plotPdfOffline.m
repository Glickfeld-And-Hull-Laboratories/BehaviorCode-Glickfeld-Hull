function [ output ] = plotPdfOffline(subjNum, dateString, blockNum, doSavePdf)
%plotMatHist: given a .mat file of a subject (subjNum) and date (dateStrng in the form "yymmdd"), this function
%creates the daily plot .pdf corresponding to it

if nargin < 3, blockNum = 'max'; end
if nargin < 4, doSavePdf = false; end

cs = holdanddetect_constants;

%dataName= strcat('data-', subjNum, '-', dateStrng, '.mat');
dataName = sprintf('data-i%03d-%s.mat', subjNum, dateString);
fileName= fullfile(cs.dataPath, dataName);

ds=mwLoadData(fileName, blockNum);

figure;
plotOnlineHist(ds, ds);

bH = findobj(gcf, 'Type', 'uicontrol');
assert(~isempty(bH));
delete(bH)

if doSavePdf
    outName = '~/Desktop/plotPdfOffline.pdf';
    exportfig_print(gcf, outName, ...
                 'FileFormat', 'pdf', ...
                 'Size', [12 12], ...
                 'PrintUI', false);
end
