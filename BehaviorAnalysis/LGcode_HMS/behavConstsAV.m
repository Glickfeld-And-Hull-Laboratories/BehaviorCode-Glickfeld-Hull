function rc = behavConstsAV


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}
        rc.pathStr = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\Behavior\Data';
        rc.dataPat = 'data-i%03d-%s.mat';
        
        if tUsername(1:5) == 'linds'
            rc.name = 'linds';
            rootDir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\lindsey';
            rc.indexFilename = fullfile(rootDir, 'Analysis\Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.fitOutputFilename = fullfile(rootDir, 'Analysis\Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDir, 'Analysis\Behavior\output\summary');
            rc.fitOutputPdfDir = fullfile(rootDir, 'Analysis\Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDir, 'Analysis\Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDir,'Analysis\Behavior\EyeTracking');
%             rc.eyeInputDir = '\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\';
            rc.lindseyAnalysis = fullfile(rootDir,'Analysis\2P');
            rc.data = fullfile(rootDir,'Data\2P_images');
            rc.ashleyAnalysis = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley\Analysis';
            rc.ashleyData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley\data';
        elseif tUsername(1:5) == 'ashle'
            rc.name = 'ashle';
            rootDirCrash = '\\CRASH.dhe.duke.edu\data\home\ashley';
%             rootDirCrashAnalysis = fullfile(rootDirCrash,'Analysis');
            rootDirIsilon = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley';
            rc.behavData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\Behavior\Data';
            rc.ashley = rootDirIsilon;
            rc.indexFilename = fullfile(rootDirIsilon,'Analysis', '\Behavior\experimentIndexes\av-subj-days-lg.xls');
            rc.indexFilename_audCtrl = fullfile(rootDirIsilon,'Analysis', '\Behavior\experimentIndexes\v-subj-days-lg.xls');
            rc.indexFilename_fsavImg = fullfile(rootDirIsilon,'Analysis', '\Behavior\experimentIndexes\fsav-img-subj-days-lg.xls');
            rc.eaMouseIndexFilename = fullfile(rootDirIsilon,'Analysis','Behavior','experimentIndexes','FSAV_eaMouseData.xls');
            rc.fitOutputFilename = fullfile(rootDirIsilon,'Analysis', '\Behavior\experimentIndexes\av-subj-fits-lg.xls');
            rc.fitOutputSummary = fullfile(rootDirIsilon,'Analysis', '\Behavior\output\analysis\summary');
            rc.fitOutputPdfDir = fullfile(rootDirIsilon,'Analysis', '\Behavior\output\pdfFits');
            rc.fitOutputMatDir = fullfile(rootDirIsilon,'Analysis', '\Behavior\output\fitMatStats');
            rc.eyeOutputDir = fullfile(rootDirIsilon,'Analysis', 'FSAV Summaries','pupil');
            rc.eyeInputDir = rootDirIsilon;
            rc.caOutputDir = fullfile(rootDirIsilon,'Analysis', 'FSAV Summaries');
            rc.ashleyAnalysis = fullfile(rootDirIsilon,'Analysis');
            rc.ashleyData = fullfile(rootDirIsilon,'data');
%             rc.isilon = rootDirIsilon;
        elseif tUsername(1:5) == 'carol'
            rc.name = 'carolyn';
            rootDir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff';
%             rc.behavData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\Behavior\Data';
            rc.behavData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\Behavior\Data';
            rc.ashleyAnalysis = fullfile(rootDir,'home\ashley\Analysis');     
            rc.ashleyData = fullfile(rootDir,'home\ashley\data');
            rc.carolynAnalysis = fullfile(rootDir,'home\carolyn\Analysis');
            rc.carolynData = fullfile(rootDir,'home\carolyn\Data');
        end
    case{'nb-hubel'}
        if tUsername(1:5) == 'cc735'
            rc.name = 'celine';
            rootDir = 'Z:\home\ACh';
            rc.celineAnalysis  = fullfile(rootDir,'Analysis\2p_analysis');
            rc.achAnalysis = fullfile(rootDir,'Analysis\2p_analysis');
            rc.achData = 'Z:\home\ACh\Data';
            
            
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

