function fitXd = readFitOutputFile
% helper function
% Read XLS/frm file, sanitize fields, check for errors
%
% histed 121226

% load file
rc = behavConstsHADC8;
fitXd = frm_xls2frm(rc.fitOutputFilename, [], rc.fitOutputTextCols);

% error checkinghn
if all(isnumeric(fitXd.DateStr))
    fitXd.DateStr = cellstr(int2str(fitXd.DateStr(:)))';
end
if isnumeric(fitXd.DateTimeStarted)
    assert(all(isnan(fitXd.DateTimeStarted)), ...
        'if numeric must be all empty - timestamp is str');
    fitXd.DateTimeStarted = cell([1 fitXd.nRows]);
end
