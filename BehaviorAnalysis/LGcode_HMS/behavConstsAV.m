function rc = behavConstsAV


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        
        if tUsername(1:5) == 'linds'
            rc.name = 'lindsey';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis';
            rc.indexFilename = fullfile(rootDir, 'Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'Behavior\output\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir,'Behavior\EyeTracking');
            rc.eyeInputDir = fullfile(rootDir,'Behavior\EyeTracking');
            rc.caOutputDir = fullfile(rootDir,'2P');
            rc.ashleyAnalysis = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashley';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis';
            rc.indexFilename = fullfile(rootDir, '\Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, '\Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, '\Behavior\output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, '\Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, '\Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir);
            rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.caOutputDir = fullfile(rootDir);
            rc.ashleyAnalysis = fullfile(rootDir);
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

