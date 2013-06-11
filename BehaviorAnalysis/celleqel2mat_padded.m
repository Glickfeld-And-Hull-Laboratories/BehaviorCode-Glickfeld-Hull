function matOut = celleqel2mat_padded (cellIn, padVal, castOutputType)
%CELLEQEL2MAT_PADDED (ps-utils): Convert cell mat of equal size els to matrix
%   MATOUT = CELLEQEL2MAT_PADDED (CELLIN, PADVAL)
%   cellIn - Each element must be the same size as every other or empty
%   padval defaults to NaN
%   castOutputType: type of output matrix
%       'double' {default}; 'int32', 'int64' 
%
%   Note: the dimension to concatenate is given by the shape of the
%   cell array.  A cell that is size [2,1] will be concatenated along
%   the first dim, and one that is [1,2] along the second.
%
%   This is MUCH faster than CELL2MAT_PADDED when each entry is either the
%   same size or empty.  It's still not terribly fast.
%
%  MH - http://github.com/histed/tools-mh

if nargin < 2, padVal = NaN; end
if nargin < 3, castOutputType = 'double'; end

if isempty(cellIn), matOut = []; return; end

emptyIx = cellfun(@isempty, cellIn);
if all(emptyIx)
  matOut = repmat(padVal, size(cellIn));
  return
end

firstNonN = find(~emptyIx, 1, 'first');
tSize = size(cellIn{firstNonN});
padMat = repmat(padVal, tSize);

[cellIn{emptyIx}] = deal(padMat);

w = warning('off', 'MATLAB:nonIntegerTruncatedInConversionToChar');

castFH = str2func(castOutputType);
matOut = cellfun(castFH, cellIn);

warning(w);

outSize = tSize .* size(cellIn);
matOut = reshape(matOut, outSize);

