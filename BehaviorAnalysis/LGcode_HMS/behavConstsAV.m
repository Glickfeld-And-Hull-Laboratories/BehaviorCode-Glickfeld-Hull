function rc = behavConstsAV


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        
        if tUsername(1:5) == 'linds'
            rc.name = 'linds';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis';
            rc.indexFilename = fullfile(rootDir, 'Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'Behavior\output\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir,'Behavior\EyeTracking');
%             rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.caOutputDir = fullfile(rootDir,'2P');
            rc.ashleyAnalysis = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley_temp';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashle';
            rootDirCrash = '\\CRASH.dhe.duke.edu\data\home\ashley';
            rootDirCrashAnalysis = fullfile(rootDirCrash,'Analysis');
            rootDirIsilon = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley_temp';
            rc.behavData = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
            rc.ashley = '\\CRASH.dhe.duke.edu\data\home\ashley';
            rc.indexFilename = fullfile(rootDirCrashAnalysis, '\Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.indexFilename_audCtrl = fullfile(rootDirCrashAnalysis, '\Behavior\experimentIndexes\v-subj-days-lg.xls');
            rc.indexFilename_fsavImg = fullfile(rootDirCrashAnalysis, '\Behavior\experimentIndexes\fsav-img-subj-days-lg.xls');
            rc.eaMouseIndexFilename = fullfile(rootDirCrashAnalysis,'Behavior','experimentIndexes','FSAV_eaMouseData.xls');
            rc.fitOutputFilename = fullfile(rootDirCrashAnalysis, '\Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDirCrashAnalysis, '\Behavior\output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDirCrashAnalysis, '\Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDirCrashAnalysis, '\Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDirCrashAnalysis, 'FSAV Summaries','pupil');
            rc.eyeInputDir = rootDirCrashAnalysis;
            rc.caOutputDir = fullfile(rootDirCrashAnalysis, 'FSAV Summaries');
            rc.ashleyAnalysis = fullfile(rootDirCrashAnalysis);
            rc.ashleyData = fullfile(rootDirCrash,'data');
        elseif tUsername(1:5) == 'xiaor'
            rc.name = 'xiaor';
            rootDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis';
            rc.behavData = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
            rc.indexFilename = fullfile(rootDir, '\Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.indexFilename_audCtrl = fullfile(rootDir, '\Behavior\experimentIndexes\v-subj-days-lg.xls');
            rc.indexFilename_fsavImg = fullfile(rootDir, '\Behavior\experimentIndexes\fsav-img-subj-days-lg.xls');
            rc.eaMouseIndexFilename = fullfile(rootDir,'Behavior','experimentIndexes','FSAV_eaMouseData.xls');
            rc.fitOutputFilename = fullfile(rootDir, '\Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, '\Behavior\output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, '\Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, '\Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir, 'FSAV Summaries','pupil');
            rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.caOutputDir = fullfile(rootDir, 'FSAV Summaries');
            rc.ashleyAnalysis = fullfile(rootDir);     
            rc.ashleyData = '\\CRASH.dhe.duke.edu\data\home\ashley\data';
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

