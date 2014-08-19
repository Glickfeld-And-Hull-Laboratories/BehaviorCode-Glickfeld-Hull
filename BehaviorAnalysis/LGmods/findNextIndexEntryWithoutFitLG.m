function [fitXd fitRowN indexRowN] = findNextIndexEntryWithoutFit; %(varargin)
%  read subj xls, find a record
% 
% histed 121003

%% arg processing
%userDefs = {};
%uo = stropt2struct(stropt_defaults(userDefs, varargin));
rc = behavConstsHADC8;

if ~exist(rc.fitOutputFilename, 'file')
    error('Output fit XLS file not found, create it with no non-header rows to start anew');
end

xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
fitXd = frm_xls2frm(rc.fitOutputFilename, [], rc.fitOutputTextCols);
%assert(length(fitXd.DateTimeStarted) > 2, 'bug - datetimestarted field');
 
fitRowN = [];  indexRowN = []; % this is returned if no data found
for iD = 1:xd.nRows
    if xd.DoFit(iD) == 0
        % force skip
        continue
    end
    foundFitRowN = frm_findrownum(fitXd, ...
        { 'Subject', xd.Subject(iD), 'DateStr', xd.DateStr{iD}, 'DataBlock', xd.DataBlock{iD} }, ...
        'IgnoreMissing', true);
    if isnan(foundFitRowN)
        % found an index field with no fit
        indexRowN = iD;
        fitRowN = fitXd.nRows+1;
        
        fitXd = frm_addrow(fitXd);
        
        fitXd.Subject(fitRowN) = xd.Subject(iD);
        fitXd.DateStr{fitRowN} = xd.DateStr{iD};
        fitXd.DataBlock{fitRowN} = xd.DataBlock{iD};
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
