function rc = behavConstsAttn


tHostname = lower(hostname);

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\CRASH.dhe.duke.edu\data\home\andrew\HMSdata';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\HMS\Behavior';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes\attn-subj-days-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output\analysis\summary');
        rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
        rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
    case {'zbook'}
        rc.pathStr = '/Users/lindsey/Desktop/Data';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '/Users/lindsey/Documents/Behavior_data';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/attn-subj-days-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output/analysis/summary');
        rc.fitOutputPdfDir = fullfile(rootDir, 'output/pdfFits');
        rc.fitOutputMatDir = fullfile(rootDir, 'output/fitMatStats');
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

