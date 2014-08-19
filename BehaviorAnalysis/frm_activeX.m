function [rawOut, typeMat] = frm_activeX(xlsFileName, sheetNum)
% FRM_ACTIVEX: Replaces frm_xlsreadpoi for WIN64 compatibility
%   [rawOut, typeMat] = FRM_ACTIVEX(xlsFileName)
%
%   rawOut is a cell array with the contents of each cell as a matlab object
%   sheetNum is the sheet number, 1-origin (first sheet numbered 1 not zero)
%   typeMat is a numeric matrix of the same size (see xlsconstantsMH for
%       mapping from values to text).  Empty cells get CELL_TYPE_BLANK
%
%   See also FRM_*
%   LG - 140819

if nargin < 2 || isempty(sheetNum), sheetNum = 1; end

if ~exist(xlsFileName, 'file')
    error('XLS file not found: %s', xlsFileName);
end

Excel = actxserver('Excel.Application');
Workbook = invoke(Excel.Workbooks,'Open', xlsFileName);
sheet = get(Excel.Worksheets,'Item', sheetNum);
invoke(sheet, 'Activate');

%find used cells by row
rowLim = 50;
rowLoc = [xlcolumn(1) '1:' xlcolumn(rowLim) int2str(1)];
ResultValue  = get(Excel.Activesheet, 'Range',rowLoc);
ResultValue.Select;
rowV = Excel.Selection.Value;
charmat = zeros(size(rowV));
for i = 1:rowLim
    charmat(1,i) = ischar(cell2mat(rowV(i)));
end
minR = min(find(charmat==0),[],2)-1;

%find used cells by column
colLim = 10000;
colLoc = [xlcolumn(1) '2:' xlcolumn(1) int2str(colLim)];
ResultValue  = get(Excel.Activesheet, 'Range',colLoc);
ResultValue.Select;
colV = Excel.Selection.Value;
nanmat = zeros(size(colV));
for i = 1:colLim-1
    nanmat(i,1) = isnan(cell2mat(colV(i)));
end
minC = min(find(nanmat==1),[],1);

if or(isempty(minR),isempty(minC))
    error('no room left in spreadsheet- change limits or make new one');
end

%pull out filled values
Location = [xlcolumn(1) '1:' xlcolumn(minR) int2str(minC)];
ResultValue  = get(Excel.Activesheet, 'Range',Location);
ResultValue.Select;
bV = Excel.Selection.Value;

%bookkeeping
siz = size(bV);
typeMat = zeros(siz);
for y = 1:siz(1)
    for x = 1:siz(2)
    typeMat(y,x) = ischar(cell2mat(bV(y,x)));
    end
end
rawOut = bV;

Excel.Quit
Excel.delete
end