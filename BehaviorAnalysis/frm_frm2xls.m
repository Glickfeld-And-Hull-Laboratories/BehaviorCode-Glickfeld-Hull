function frm_frm2xls(dsT, xlsFileName, sheetNumOrStr)
%FRM_FRM2XLS: write data frame from FRM_XLS2FRM back to Excel file
%
%   FRM_FRM2XLS(dsT, xlsFileName, sheetNumOrStr)
%
%   Must call FRM_JAVASETUP first to set up java path w/ POI jar files
%
%   See also FRM_*
%
%  MH - http://github.com/histed/tools-mh


if nargin < 3, sheetNumOrStr = 1; end

%% massage colnames for output
fNames = fieldnames(dsT);
removeIx = ismember(fNames, {'colNames', 'nRows', 'nCols'});
fNames = fNames(~removeIx);
tColNames = fNames;
tColNames = regexprep(tColNames, '([A-Z])', ' $1'); % add spaces before each capital
tColNames = regexprep(tColNames, '([0-9]*)', ' $1');  % add spaces before runs of numbers
tColNames = regexprep(tColNames, '^ ', '');  % remove any leading spaces; easier to do this than use lookbehind

%% create the raw cell matrix
nCols = length(fNames);
nRows = length(dsT.(fNames{1}));
rawCell = cell(nRows+1, nCols);
rawCell(1,:) = tColNames;

for iC = 1:nCols
    tVals = dsT.(fNames{iC});
    if isnumeric(tVals)
        tValsO = mat2cell_singleton(tVals);
        [tValsO{isnan(tVals)}] = deal([]);
        tVals = tValsO;
    end   % else we keep as a cell vector
    if length(tVals) ~= nRows
        error('malformed data struct - each field must have same number of rows');
    end       
    rawCell(2:end,iC) = tVals(:);
end

%% write
 if ispc
    xlswrite(xlsFileName, rawCell, sheetNumOrStr);
 elseif ismac
    xlwrite(xlsFileName, rawCell, sheetNumOrStr);
 else
     error('not pc or mac- need to figure out how to write xls')
 end

