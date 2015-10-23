function rc = behavConstsAV


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        
        if tUsername(1:5) == 'linds'
            rc.name = 'lindsey';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Behavior';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'output\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir,'EyeTracking');
            rc.eyeInputDir = fullfile(rootDir,'EyeTracking');
            rc.caOutputDir = fullfile(rootDir,'2P');
            rc.ashleyAnalysis = '\\CRASH.dhe.duke.edu\data\home\ashley\analysis';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashley';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\Behavior';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir,'Summaries\eye tracking');
            rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
       elseif tUsername(1:5) == 'aiwei'
   
            rc.name = 'aiwei';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\aiwei\Analysis\Behavior_data';
            rc.indexFilename = fullfile(rootDir, 'experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\av-subj-fits-lg.xls');
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

