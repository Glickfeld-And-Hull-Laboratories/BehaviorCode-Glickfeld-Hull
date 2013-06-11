function frm_xlswritepoi(xlsFileName, rawCell, sheetNumOrStr, doMatInCellAsText)
%FRM_XLSWRITEPOI: use Apache POI Java libraries to write excel files w/ type info etc.
%
%   Note - empty elements of rawCell are untouched if excel file exists.   
%   
%   Must call xlsjavasetupMH.m first to set up java path w/ POI jar files
%
%   See also FRM_*
%
%  MH - http://github.com/histed/tools-mh

% notes: To include formulas, they will likely have to be stored as strings and 
% we will have to pass in the type mat from frm_xlsreadpoi.m

if nargin < 3, sheetNumOrStr = 1; end
if nargin < 4, doMatInCellAsText = true; end

import org.apache.poi.ss.usermodel.*;

try
    b0 = Cell.CELL_TYPE_NUMERIC;
catch
    error('Java classes not found - did you run xlsjavasetupMH?');
end

%% set up
xc = frm_constants;
[xPath, xName, xExt] = fileparts(xlsFileName);
if isempty(xExt), xExt = '.xlsx'; end
xlsFileName = fullfile(xPath, [xName xExt]);

if any(xlsFileName == '~')
    error('Tilde expansion not supported by Java POI library: use USER var or specify full dir');
end

%% read file if existing, else create
if exist(xlsFileName, 'file')
    jFH = java.io.FileInputStream(xlsFileName);
    wb = WorkbookFactory.create(jFH);
else
    % create it - always xlsx
    %wb = org.apache.poi.xssf.usermodel.XSSFWorkbook();
    wb = org.apache.poi.hssf.usermodel.HSSFWorkbook(); % xls
end

% debug
%fileOut = java.io.FileOutputStream(xlsFileName);
%wb.write(fileOut);
%fileOut.close();

%% find sheet if existing, else create it
tS = [];
if isnumeric(sheetNumOrStr)
    nSheets = wb.getNumberOfSheets();
    if sheetNumOrStr <= nSheets
        tS = wb.getSheetAt(sheetNumOrStr-1);
    end
elseif ischar(sheetNumOrStr)
    tS = wb.getSheet(sheetNumOrStr);  % returns empty if not found
end
if isempty(tS)
    if isnumeric(sheetNumOrStr)
        sheetName = sprintf('Sheet%02d', sheetNumOrStr);
    else
        sheetName = sheetNumOrStr;
    end
    tS = wb.createSheet(sheetName);
end

%% compute types
%  we only use string and numeric for now (see above for formulas)
[nRows nCols] = size(rawCell);

charIx = cellfun(@ischar, rawCell);
numericIx = cellfun(@isnumeric, rawCell);
emptyIx = cellfun(@isempty, rawCell);

%% iterate over rows and cols and write data
for iR = 1:nRows
    tRow = tS.getRow(iR-1);
    if isempty(tRow)
        tRow = tS.createRow(iR-1);
    end
    
    for iC = 1:nCols
        tVal = rawCell{iR,iC};
        if isempty(tVal)
            continue  % no value, so don't make any Java calls this iteration
        else
            tCell = tRow.getCell(iC-1);
            if isempty(tCell)
                tCell = tRow.createCell(iC-1);
            end
            
            if ischar(tVal)
                tCell.setCellType(Cell.CELL_TYPE_STRING);
                tCell.setCellValue(tVal);
            elseif isnumeric(tVal) && length(tVal) == 1
                tCell.setCellType(Cell.CELL_TYPE_NUMERIC);
                tCell.setCellValue(tVal);
            elseif isnumeric(tVal) && length(tVal) > 1
                if doMatInCellAsText
                    tCell.setCellType(Cell.CELL_TYPE_STRING)
                    tCell.setCellValue(mat2str(tVal));
                else
                    error('Matrix found in cell, and doMatInCellAsText is false');
                end
            else
                error('unknown type');
            end
        end
    end
end

%% write file to disk and close
fileOut = java.io.FileOutputStream(xlsFileName);
wb.write(fileOut);
fileOut.close();
