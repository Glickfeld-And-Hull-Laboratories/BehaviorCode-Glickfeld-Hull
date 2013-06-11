function xd = readIndexFile
% helper function
% Read XLS/frm file, sanitize fields, check for errors
%
% histed 121226

% load file
rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

chkUniqueId(xd.UniqueId);

% do explicit casting for cases where xls may have different types
if all(isnumeric(xd.DateStr))
    xd.DateStr = cellstr(int2str(xd.DateStr(:)))';
end

% should do strtrim in here too on all the index fields
