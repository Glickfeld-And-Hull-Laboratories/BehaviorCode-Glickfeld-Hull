function dsT = xlsread_textornum(xlsFileName, sheet)
%XLSREADFRAMEMH: read xls with header row; convert to struct by col title
% 
%   ds = xlsread_textornum(xlsFileName, sheet)
%
%   ds is a structure.  The first row is treated as a header line and the
%   cell contents are used as field names in the structure. 
%   Each field is a vector from the column contents.  If there is mixed
%   numeric and string data, a raw cell vector is returned.  If only
%   numeric data we return a numeric vector, else a text vector.
%
%   Uses MATLAB XLSREAD function, not Apache POI Java library
%
%  MH - http://github.com/histed/tools-mh

% histed 120530: created

if nargin < 2, sheet = ''; end

warning('off', 'MATLAB:xlsread:Mode');
[~,~,raw] = xlsread(xlsFileName, sheet);

% remove header line
dsT.colNames = raw(1,:);
raw = raw(2:end,:);  
% convert to a struct
dsT.nCols = length(dsT.colNames);
dsT.nRows = size(raw,1);

isSomeNum = any(cellfun(@isnumeric,raw), 1);
isSomeText = any(cellfun(@ischar,raw), 1);

for iC=1:dsT.nCols
    tFN = dsT.colNames{iC};
    % sanitize fname
    if isnan(tFN)
        tFN = sprintf('Column%02d', iC);
    else
        tFN = strrep(tFN, sprintf('\n'), '');
        tFN = regexprep(tFN, '[-\(\)]', '_');  % misc punct w/ underscores
        tFN = genvarname(tFN);
        tFN = regexprep(tFN, '_$', '');  % misc punct w/ underscores
    end
    
    if isSomeNum(iC) && isSomeText(iC)
        tV = raw(:,iC);
    elseif isSomeNum(iC)
        tV = celleqel2mat_padded(raw(:,iC), NaN);
    else % text
        tV = cellstr(raw(:,iC));
    end
        
    dsT.(tFN) = tV(:)';
end
