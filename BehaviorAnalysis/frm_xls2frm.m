function dsT = frm_xls2frm(xlsFileName, sheet, textFieldNames)
%FRM_XLS2FRM: read xls with header row; convert to struct by col title
% 
%   ds = FRM_XLS2FRM(xlsFileName, sheet, textFieldNames)
%
%   ds is a structure.  The first row is treated as a header line and the
%   cell contents are used as field names in the structure. 
%   Each field is a vector from the column contents.  If there is mixed
%   numeric and string data, a raw cell vector is returned.  If only
%   numeric data we return a numeric vector, else a text vector.
%
%   if textFieldNames is specified, convert the columns with these names to text even if found 
%      to be all numeric data
%
%   Uses Apache POI Java library for Excel file processing
%
%   See also: FRM_*
% 
%  MH - http://github.com/histed/tools-mh

% histed 120530: created

if nargin < 2, sheet = 1; end
if nargin < 3, textFieldNames = []; end

xc = frm_constants;

%[raw, typeMat] = frm_xlsreadpoi(xlsFileName, sheet); %old java code from MH 
[raw, typeMat] = frm_activeX(xlsFileName, sheet); %new code- LG 140819

% remove header line
dsT.colNames = raw(1,:);
raw = raw(2:end,:);  
% convert to a struct
dsT.nCols = length(dsT.colNames);
dsT.nRows = size(raw,1);

isSomeNum = any(typeMat(2:end,:) == xc.typeNums.NUMERIC,1);
isSomeText = any(typeMat(2:end,:) == xc.typeNums.STRING,1);
emptyIx = cellfun(@isempty, raw(:,:));

removeColsIx = false(size(dsT.colNames));
for iC=1:dsT.nCols
    tFN = dsT.colNames{iC};
    % sanitize fname
    if isempty(tFN)
        assert(all(emptyIx(:,iC)), 'empty column name but data in rows below');
        removeColsIx(iC) = true;
        continue;
    elseif isnan(tFN)
        tFN = sprintf('Column%02d', iC);
    else
        tFN = strrep(tFN, sprintf('\n'), '');
        tFN = regexprep(tFN, '[-\(\)]', '_');  % misc punct w/ underscores
        tFN = genvarname(tFN);
        tFN = regexprep(tFN, '_$', '');  % misc punct w/ underscores
    end
    %if strcmp(tFN, 'DateTimeStarted'), keyboard, end
    %if strcmp(tFN, 'b2OneRampExtraConstantLengthMs'), keyboard, end

    if isSomeNum(iC) && isSomeText(iC)
        % hack for now - convert all numeric elements to text MH 120801
        % MH 130125: This mainly works, we may want to do a better job of allowing caller 
        %  to do explicit casting (ie. the textFieldNames input)

        tV = raw(:,iC);
        numIx = typeMat(2:end,iC) == xc.typeNums.NUMERIC;
        desIx = numIx & ~emptyIx(:,iC);
        if sum(desIx)>0
            tV(desIx) = cellstr(num2str(cat(1,tV{desIx})));
            %warning('converting partial num/text field to all text: right decision?');
        end
        
        % convert empty numeric matrices to string
        desIx = cellfun(@isempty, tV);
        [tV{desIx}] = deal('');
    elseif isSomeNum(iC) && ~isSomeText(iC)
        tV = celleqel2mat_padded(raw(:,iC), NaN);
    elseif ~isSomeNum(iC) && ~isSomeText(iC)
        % empty
        tV = nan(size(raw,1),1);
    else % ~isSomeNum && isSomeText
        tV = raw(:,iC);
        eIx = emptyIx(:,iC);
        [tV{eIx}] = deal('');    % convert double empties to char 
        %tV = cellstr(tV);  % this trims blanks, so just leave as straight cell
    end
        
    dsT.(tFN) = tV(:)';
end

% explicit text conversion
if ~isempty(textFieldNames)
    nF = length(textFieldNames);
    for iF = 1:nF
        tFN = textFieldNames{iF};
        if ~isfield(dsT, tFN)
            error('requested textFieldName %s not found in data', tFN);
        end
        if isnumeric(dsT.(tFN))
            % convert
            dsT.(tFN) = cellstr(num2str(dsT.(tFN)'))';
        else
            assert(iscell(dsT.(tFN)), 'fields must be numeric or cell');
        end
    end
end

        
            

dsT.colNames = dsT.colNames(~removeColsIx);
dsT.nCols = length(dsT.colNames);
