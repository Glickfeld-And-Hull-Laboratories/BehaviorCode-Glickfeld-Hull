function rc = behavConstsHADC8


tHostname = lower(hostname);

switch tHostname
    case {'zbook'}
        rc.pathStr = '/Users/lindsey/Desktop/Data';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '/Users/lindsey/Documents/Behavior_data';
        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-days-lg.xls');
        rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-fits-lg.xls');
        rc.fitOutputSummary = fullfile(rootDir, 'output/analysis/summary');

    case {'mambo'}
        rc.pathStr = '/Users/histed/data/mus-behavior/SnowLeopard/Data';
        rc.dataPat = 'data-i%03d-%s.mat';

        rootDir = '/Users/histed/analysis/mus-project/ramptrain-analysis1';

        rc.indexFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-days.xls');
        rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes/chr2-subj-fits.xls');

end

rc.fitOutputTextCols = { 'DateStr', 'DataBlock', 'DateTimeStarted', 'PdfFigFilename', 'MatFilename' };
rc.indexTextCols = { 'DateStr', 'DataBlock', 'TrialRangeToUse', 'Notes' };

rc.fitOutputPdfDir = fullfile(rootDir, 'output/pdfFits');
rc.fitOutputMatDir = fullfile(rootDir, 'output/fitMatStats');
rc.fhWeibull = @(p,xs) (p(4) + (p(3)-p(4)) .* (1 - exp( -(xs./p(1)) .^ p(2))));


%%%%%%%%%%% simple simple functions

rc.computeFName = @subComputeFName;

function outFName = subComputeFName(subjNum, dateStr)
    outFName = fullfile(rc.pathStr, sprintf(rc.dataPat, subjNum, deblank(dateStr)));
end



end

