function [rawOut, typeMat] = frm_xlsreadpoi(xlsFileName, sheetNum)
%FRM_XLSREADPOI: use Apache POI Java libraries to read excel files w/ type info etc.
%   [rawOut, typeMat] = FRM_XLSREADPOI(xlsFileName)
%
%   rawOut is a cell array with the contents of each cell as a matlab object
%   sheetNum is the sheet number, 1-origin (first sheet numbered 1 not zero)
%   typeMat is a numeric matrix of the same size (see xlsconstantsMH for
%       mapping from values to text).  Empty cells get CELL_TYPE_BLANK
%
%   Must call frm_javasetup.m first to set up java path w/ POI jar files
%
%   See also FRM_*
%
%  MH - http://github.com/histed/tools-mh

% v2, 120801: move iteration to java code to speed it up.  (each java call in
% matlab is very slow)

if nargin < 2 || isempty(sheetNum), sheetNum = 1; end

if ~exist(xlsFileName, 'file')
    error('XLS file not found: %s', xlsFileName);
end

tObj = frm_jcReadXls(xlsFileName,sheetNum);  % java
bV = tObj.getCellContents();

% unpack outputs
numMat = bV(1);
strCell = cell(bV(2));
typeMat = double(bV(3));

% make a single cell array
isNum = ~isnan(numMat);
isStr = ~cellfun(@isempty, strCell);
assert(~any(any(isNum & isStr)));
rawOut = strCell;
numCell = mat2cell_singleton(numMat);
rawOut(isNum) = numCell(isNum);

%typeMat = celleqel2mat_padded(typeC, Cell.CELL_TYPE_BLANK);
