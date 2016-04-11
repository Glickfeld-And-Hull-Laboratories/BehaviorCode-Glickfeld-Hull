function rc = behavConstsHADC8


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'zbook'}
        rc.pathStr = '/Users/lindsey/Desktop/Data';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '/Users/lindsey/Documents/Behavior_data';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-days-lg.xls');
        rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-fits-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output/analysis/summary');
        rc.fitOutputPdfDir = fullfile(rootDir, 'output/pdfFits');
        rc.fitOutputMatDir = fullfile(rootDir, 'output/fitMatStats');
    
    case {'jeffbook pro'}
        rc.pathStr = '/Users/jeffreysims/Documents/Behavior Data/Data';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '/Users/jeffreysims/Documents/Behavior Data';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/subj-days-lg.xls');
        rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes/subj-fits-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output/analysis/summary');
        rc.fitOutputPdfDir = fullfile(rootDir, 'output/pdfFits');
        rc.fitOutputMatDir = fullfile(rootDir, 'output/fitMatStats');

    case {'lgair'}
        rc.pathStr = 'andrew/Behavior/Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        rootDir = '/Users/llglick/Desktop';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/subj-days-lg.xls');
        rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes/subj-fits-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output/analysis/summary');
        rc.fitOutputPdfDir = fullfile(rootDir, 'output/pdfFits');
        rc.fitOutputMatDir = fullfile(rootDir, 'output/fitMatStats');
        
    case {'nuke'}
        if tUsername(1:6) == 'lindse'
            rc.pathStr = 'Z:\home\andrew\Behavior\Data';
            rc.dataPat = 'data-i%03d-%s.mat';
            rootDir = 'Z:\home\lindsey\Analysis\Behavior';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'output\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
        elseif tUsername(1:6) == 'miaomi'
            rc.pathStr = 'Z:\andrew\Behavior\Data';
            rc.dataPat = 'data-i%03d-%s.mat';
            rootDir = 'Z:\miaomiao\behavior\Behavior';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\subj-days-mj.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\subj-fits-mj.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
        elseif tUsername(1:6) == 'ashley'
            rc.pathStr = 'Y:\home\andrew\Behavior\Data';
            rc.dataPat = 'data-i%03d-%s.mat';
            rootDir = 'Y:\home\ashley\Analysis\Behavior';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
        end
        
end

rc.fitOutputTextCols = { 'DateStr', 'DataBlock', 'DateTimeStarted', 'PdfFigFilename', 'MatFilename' };
rc.indexTextCols = { 'DateStr', 'DataBlock', 'TrialRangeToUse', 'Notes' };

rc.fhWeibull = @(p,xs) (p(4) + (p(3)-p(4)) .* (1 - exp( -(xs./p(1)) .^ p(2))));


%%%%%%%%%%% simple simple functions

rc.computeFName = @subComputeFName;

function outFName = subComputeFName(subjNum, dateStr)
    outFName = fullfile(rc.pathStr, sprintf(rc.dataPat, subjNum, deblank(dateStr)));
end



end

