function rc = behavConstsDART


tHostname = lower(hostname);
[s,tUsername] = dos('ECHO %USERNAME%');

switch tHostname
    case {'nuke'}   
        if tUsername(1:5) == 'linds'
            rc.name = 'linds';
            rootDir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\lindsey';
            rc.analysis = fullfile(rootDir,'Analysis\2P');
            rc.data = fullfile(rootDir,'Data\2P_images');
            rc.ashleyAnalysis = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley\Analysis';
            rc.ashleyData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\ashley\data';
            rc.tammyAnalysis = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\tammy\Analysis';
            rc.tammyData = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\all_staff\home\tammy\data';
        elseif tUsername(1:5) == 'tammy'
            rc.name = 'tammy';
            rootDir = 'Z:\home\ACh';
            rc.achAnalysis = fullfile(rootDir,'Analysis','2p_analysis');
            rc.achData = 'Z:\home\ACh\Data';
            rc.analysis = fullfile(rootDir,'Analysis');
            rc.data = 'Z:\home\ACh\Data';

        elseif tUsername(1:5) == 'celine'
            rc.name = 'celine';
            rootDir = 'Z:\home\ACh';
            rc.achAnalysis = fullfile(rootDir,'Analysis');
            rc.achData = 'Z:\home\ACh\Data';
            rc.analysis = fullfile(rootDir,'Analysis');
            rc.data = 'Z:\home\ACh\Data';
           
        end
    case{'nb-hubel'}
        if tUsername(1:5) == 'cc735'
            rc.name = 'celine';
            rootDir = 'Z:\home\ACh';
            rc.achAnalysis = fullfile(rootDir,'Analysis','2p_analysis');
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

