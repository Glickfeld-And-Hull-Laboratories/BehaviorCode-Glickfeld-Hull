function rowN = frm_findrownum(xd, indexCell, varargin)
%FRM_FINDROWNUM: given a set of fields serving as an index, find row
%
%   rowN = frm_findrownum(xd, indexCell, [options])
%
%   Find the unique row(s) matching the index fields.
%   To specify multiple rows make the index values vectors/cell vectors
%  
%   Options: 
%        IgnoreMissing: true: return NaN if not found; 
%                             false (default): raise error
%
%   Example: 
%      rowN = frm_findrownum(xd, { 'File', 'ab', 'Extension', 'bc' };
%
% histed 120531

userDefs = { ...
    'IgnoreMissing', false, ...
    'AllowMultiple', false, ...  
    };

uo = stropt2struct(stropt_defaults(userDefs, varargin));

debug = 0;

% firstpass check if any data
if xd.nRows == 0
    if uo.IgnoreMissing
        rowN = NaN;
        return
    else
        error('No rows in data - no rows can be matched');
    end
end

desFNames = indexCell(1:2:end);
desFValues = indexCell(2:2:end);

nFields = length(desFNames);

% convert chars to singleton cells
strIx = cellfun(@ischar, desFValues) | cellfun(@isempty, desFValues);
newC = cellfun(@(x) {{x}}, desFValues(strIx));
desFValues(strIx) = newC;

nOut = length(desFValues{1});
rowNOut = nan([1 nOut]);

%: allow empties to find unfilled cells now MH 130123
% error checking
% zeroIx = cellfun(@length, desFValues)==0;
% if any(zeroIx)
%     error('Missing value in input: %s has length zero', desFNames{find(zeroIx,1)});
% end


for iN = 1:nOut
    desIx = true([1 xd.nRows]);
    for iF=1:nFields
        tFN = desFNames{iF};
        tFV = desFValues{iF};

        % deal with arrays of inputs; take just one element
        if isnumeric(tFV) && length(tFV) > 1
            % input is a vector, take just one element
            assert(length(tFV) == nOut, 'Non matching input columns');
            tFV = tFV(iN);
        elseif iscell(tFV)
            assert(length(tFV) == nOut, 'Non matching input columns');
            tFV = tFV{iN};
        end

        assert( ( isempty(tFV) ...
            || ( (isnumeric(tFV)|isa(tFV, 'function_handle')) && length(tFV)==1) ...
            || ischar(tFV) ), ...
            'bug: only deals with singletons below');
        
        comparisonType = '';
        compFn = [];
        
        % choose comparison based on rows
        if iscell(xd.(tFN))
            comparisonType = 'cell';
        elseif isnumeric(xd.(tFN))
            comparisonType = 'equals';
        else
            error('Invalid column type');
        end
        
        % choose compfn based on input type
        if isa(tFV, 'function_handle')
            compFn = tFV;
        elseif isnan(tFV)
            compFn = @isnan;
        elseif isempty(tFV)
            compFn = @isempty;
            if isnumeric(xd.(tFN))
                error('Looking for empty in numeric array: probably want NaN instead');
            end
        elseif ischar(tFV);
            compFn = @(x) strcmp(tFV, char(x));
        else
            compFn = @(x) eq(tFV, x);  % hopefully I didn't miss a case
        end

        % do comparison
        if strcmp(comparisonType, 'cell')
            tDIx = cellfun(compFn, xd.(tFN));
        else
            tDIx = feval(compFn, xd.(tFN));
        end
        
        % check for errs
        if sum(tDIx) == 0 && ~uo.IgnoreMissing
            error('Index value matched 0 rows: Field %s, desired value %s', ...
                tFN, mat2str(tFV));
        end
        
        if debug
            nM = sum(desIx&tDIx);
            if nM <= 3
                mStr = [', rows: ' mat2str(find(desIx&tDIx))];
            else
                mStr = '';
            end
            
            fprintf(1, 'Debug: Field %s val %s, restricted from %d to %d rows (matched %d)%s\n', ...
                tFN, mat2str(tFV), sum(desIx), sum(desIx&tDIx), sum(tDIx), mStr);
        end
        desIx = desIx & tDIx;
    end
    
    nD = sum(desIx);
    if nD < 1
        if ~uo.IgnoreMissing
            outStr = evalc('columnNames = desFNames, valuesToMatch = desFValues');
            error('No index matches found: \n %s', outStr);
        else
            rowN(iN) = NaN;
            continue
        end
    elseif nD > 1
        if uo.AllowMultiple
            assert(nOut == 1, 'if AllowMultiple is true, each input field match must be length 1');
            rowN = find(desIx);
        else
            error('More than one match found and AllowMultiple is false');
        end
    elseif nD == 1
        rowN(iN) = find(desIx);
    end
end
%assert(~any(isnan(rowN)));
