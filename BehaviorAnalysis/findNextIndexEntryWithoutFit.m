function [fitXd fitRowN indexRowN] = findNextIndexEntryWithoutFit
%  read subj xls, find a record
% 
% histed 121003

%% arg processing
rc = behavConstsHADC8;

if ~exist(rc.fitOutputFilename, 'file')
    error('Output fit XLS file not found, create it with no non-header rows to start anew');
end

xd = readIndexFile;
fitXd = readFitOutputFile;
 
fitRowN = [];  indexRowN = []; % this is returned if no data found
for iD = 1:xd.nRows
    if xd.DoFit(iD) == 0
        % force skip
        continue
    end
    if isnan(xd.Subject(iD))
        error('behaviorAnalysisHADC8MH:NoMoreRowsToDo', ...
            'Stopping fits; prev rows are done and Subject field empty in row %d', iD);
    end
    foundFitRowN = frm_findrownum(fitXd, ...
        { 'Subject', xd.Subject(iD), 'DateStr', xd.DateStr{iD}, 'UniqueId', xd.UniqueId{iD} }, ...
        'IgnoreMissing', true);
    if isnan(foundFitRowN)
        % found an index field with no fit
        indexRowN = iD;
        fitRowN = fitXd.nRows+1;
        
        fitXd = frm_addrow(fitXd);
        
        fitXd.Subject(fitRowN) = xd.Subject(iD);
        fitXd.DateStr{fitRowN} = xd.DateStr{iD};
        fitXd.UniqueId{fitRowN} = xd.UniqueId{iD};
        break;
    elseif isempty(fitXd.DateTimeStarted{foundFitRowN})
        % do existing row w/ no DateTimeStarted field
        fitRowN = foundFitRowN;
        indexRowN = iD;
        break
    else
        % data is already there
        continue
    end
end
