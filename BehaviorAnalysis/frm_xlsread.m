function [rawOut, typeMat] = frm_xlsread(xlsFileName, sheetNum)
% FRM_ACTIVEX: Replaces frm_xlsreadpoi for WIN64 compatibility
%   [rawOut, typeMat] = FRM_ACTIVEX(xlsFileName)
%
%   rawOut is a cell array with the contents of each cell as a matlab object
%   sheetNum is the sheet number, 1-origin (first sheet numbered 1 not zero)
%   typeMat is a numeric matrix of the same size (see xlsconstantsMH for
%       mapping from values to text).  Empty cells get CELL_TYPE_BLANK
%
%   See also FRM_*
%   LG - 141024

if nargin < 2 || isempty(sheetNum), sheetNum = 1; end

if ~exist(xlsFileName, 'file')
    error('XLS file not found: %s', xlsFileName);
end

[a b rawOut] = xlsread(xlsFileName, sheetNum);

%bookkeeping
siz = size(rawOut);
typeMat = zeros(siz);
for y = 1:siz(1)
    for x = 1:siz(2)
    typeMat(y,x) = ischar(cell2mat(rawOut(y,x)));
    end
end

end