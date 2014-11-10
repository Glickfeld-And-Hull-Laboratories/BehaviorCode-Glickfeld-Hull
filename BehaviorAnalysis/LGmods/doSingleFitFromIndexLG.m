function [fitX1d figH fitS bootS] = doSingleFitFromIndexLG(varargin)
% read subj xls, find a record, do fits, write output to xls
%
%  DoSubjAndDate can be 'next' or { subjNumAsInt, dateStr }
% 
% histed 120602

%% arg processing
userDefs = { ...
    'DoSubjAndDate', 'next', ...
    'DataBlock', 'max', ...
    'NBootstrapReps', 1000, ...
    'LoadDataDebug', false, ...
    'DoExport', false };
    
uo = stropt2struct(stropt_defaults(userDefs, varargin));


%% get xls data and find params
rc = behavConstsHADC8;

%% select which day's data and set params
if ischar(uo.DoSubjAndDate) && strcmp(uo.DoSubjAndDate, 'next')
    [fitXd fitRowN indexRowN] = findNextIndexEntryWithoutFitLG;

    if isempty(indexRowN)
        disp('No more records to do');
        fitX1d = []; figH = []; fitS = []; bootS = [];
        return
    end

    % do the one we found
    subjNum = fitXd.Subject(fitRowN);
    dateStr = fitXd.DateStr{fitRowN};
    dataBlock = fitXd.DataBlock{fitRowN};
else
    % specified directly
    subjNum = uo.DoSubjAndDate{1};
    dateStr = uo.DoSubjAndDate{2};
    dataBlock = uo.DataBlock;
    fitXd = [];
end

%% get xls data 
fprintf(1, 'Starting subject %03d, date %s, block %s\n', ...
    subjNum, dateStr, mat2str(dataBlock));
[bs ds x1d fitX1d] = getAllBehavDataUsingIndexLG(...
    'Subject', subjNum, ...
    'DateStr', dateStr, ...
    'DataBlock', dataBlock, ...
    'LoadDataDebug', uo.LoadDataDebug, ...
    'FitXd', fitXd);  % we need to use this because the above may add a row to fitXd

%% setup fit params from index
fitParams = { ...
    'DoCorrectEarlies', x1d.DoCorrectEarlies, ...
    'DoClampAtZero', true }; % always clamp at zero, esp for 55 120526


%% do fit and plot
block2IndexList = [0 1];
figH = figure('Tag', 'doSingleFitFromIndex');%(53);
clear fitS bootS

if ds.is2AFC
    bs.doClampAtZero = 0.5;
    bs.nsides = 3;
else
    bs.doClampAtZero = true;
    bs.nsides = 1;
end
bs.ds = ds;
bs.x1d = x1d;
for iB = 1:bs.nBlock2Indices
    bs.iB = iB;
    tBN = block2IndexList(iB);
    whichTrials = bs.whichTrials(find(bs.block2V == tBN));
    if iB==1
        bs.rep = 1;
    end
    [fitS(iB) bootS(iB)] = weibullFitWithBootstrapLG(...
        'BehavDataStruct', bs, ...
        'DoPlot', true, ...
        fitParams{:}, ...
        'DoBootstrap', true, ...
        'NBootstrapReps', uo.NBootstrapReps);
end

suptitle2(sprintf('i%03d %s', subjNum, dateStr));

%% load into output mat, write
desN = frm_findrownum(fitXd, ...
    { 'Subject', subjNum, ...
    'DateStr', dateStr, ...
    'DataBlock', dataBlock });
dbStr = mat2str(dataBlock);
dbStr = regexprep(dbStr, '[\[\]]', '_');  % remove quotes, brackets, spaces commas
dbStr = regexprep(dbStr, '[''" ,]', '');  % remove quotes, brackets, spaces commas

% create missing dirs
if ~exist(rc.fitOutputMatDir, 'dir')
    mkdir(rc.fitOutputMatDir);
end
if ~exist(rc.fitOutputPdfDir, 'dir')
    mkdir(rc.fitOutputPdfDir);
end


% save all stats to a mat
outMatName = fullfile(rc.fitOutputMatDir, ...
    sprintf('subj%03d-%s-%s.mat', subjNum, dateStr, dbStr));
save(outMatName, 'fitS', 'bootS');

outName = fullfile(rc.fitOutputPdfDir, ...
    sprintf('subj%03d-%s-%s.pdf', subjNum, dateStr, dbStr));
for iB = 1:bs.nBlock2Indices
    for iS = 1;
    fitXd.(sprintf('Threshold%d', iB))(desN) = fitS(iB).thresh(iS);
    fitXd.(sprintf('Slope%d', iB))(desN) = fitS(iB).slope(iS);
    fitXd.(sprintf('Threshold%dCi95Low', iB))(desN) = fitS(iB).bootStats.ci95(1,iS);
    fitXd.(sprintf('Threshold%dCi95High', iB))(desN) = fitS(iB).bootStats.ci95(2,iS);
    fitXd.(sprintf('Threshold%dCi99Low', iB))(desN) = fitS(iB).bootStats.ci99(1,iS);
    fitXd.(sprintf('Threshold%dCi99High', iB))(desN) = fitS(iB).bootStats.ci99(2,iS);
    fitXd.(sprintf('Slope%dCi95Low', iB))(desN) = fitS(iB).bootStats.slopeCi95(1,iS);
    fitXd.(sprintf('Slope%dCi95High', iB))(desN) = fitS(iB).bootStats.slopeCi95(2,iS);
    fitXd.(sprintf('Slope%dCi99Low', iB))(desN) = fitS(iB).bootStats.slopeCi99(1,iS);
    fitXd.(sprintf('Slope%dCi99High', iB))(desN) = fitS(iB).bootStats.slopeCi99(2,iS);
    fitXd.PdfFigFilename{desN} = outName;
    fitXd.MatFilename{desN} = outMatName;
    fitXd.DateTimeStarted{desN} = datestr(now);
    end
end
fitX1d = frm_extractrow(fitXd, desN);

% write output
if uo.DoExport 
    exportfig_print(gcf, outName, ...
        'FileFormat', 'pdf', ...
        'Size', [6 4.8+(4*(bs.nBlock2Indices-1))]);
    
    frm_frm2xls(fitXd, rc.fitOutputFilename);
end
