function singleStruct = structarray2singleton(structArrIn)
%
%  Note that this uses cell2mat on each field to make it a single matrix, so
%  it is subject to the limitations of cell2mat.  
%
%  However, we first look for double matrices with empty entries; empty is
%  padded with NaN.
%
%  MH - http://github.com/histed/tools-mh

error('should use struct_vect2singleton');

anonCell = struct2cell(structArrIn);
fNames = fieldnames(structArrIn);
anonCellSize = size(anonCell);
nFields = anonCellSize(1);
structArrSize = anonCellSize(2:end);

for iF=1:nFields
    tFieldName = fNames{iF};
    
    tFieldVals = reshape(anonCell(iF,:), structArrSize); %handle: multi-dims
                                                         %(of orig struct arr)    
    isDIx = cellfun(@isfloat, tFieldVals);
    isEIx = cellfun(@isempty, tFieldVals);
    isSIx = cellfun(@isstruct, tFieldVals);
    if all(isDIx(:) | isEIx(:))
        % every cell is a double or empty, fill with NaN
        if all(isEIx(:))
            % all empty
            singleStruct.(tFieldName) = []; 
            continue;
        else
            % some el numeric, some empty
            [tFieldVals{isEIx}] = deal(NaN);
        end
        
        singleStruct.(tFieldName) = cell2mat(tFieldVals);
    elseif all(isSIx(:) | isEIx(:))
        % structure, call recursively
        
        if all(isSIx(:))
            s = structincell2structarray(tFieldVals);
            singleStruct.(tFieldName) = structarray2singleton(s);
        end
    else
        % copy directly
        singleStruct.(tFieldName) = tFieldVals;
    end

end

