function output = smoothfit1(varargin)
%SMOOTHFIT1 (tools-mh): apply SMOOTH along one dim of a matrix
%   OUTPUT = SMOOTHFIT1(..., 'dim', 1)
%
%   See also SMOOTH (curve fit toolbox)
%
%  MH - http://github.com/histed/tools-mh

% checks
if ~exist('smooth') == 2
    error('SMOOTH function not found - curve fit toolbox missing?');
end

% break out dimension arg
dimIx = strcmpi(varargin, 'dim') | strcmpi(varargin, 'dimension');
dimN = find(dimIx);

nD = length(dimN);
if nD > 1
    error('Too many dimension arguments');
end

desArgs = 2:length(varargin); % skip first mat input
mat = varargin{1};
matSize = size(mat);
assert(ndims(mat) == 2, 'Input must have only two dimensions');

% find relevant dimension to operate on
if nD == 1
    tDim = varargin{dimN+1};
    desArgs = setdiff(desArgs, [dimN dimN+1]);
    assert(matSize(tDim) > 1, 'Asked to work along a dimension of length 1');
elseif nD == 0
    if matSize(1) == 1
        tDim = 2;
    elseif matSize(2) == 1
        tDim = 1;
    else
        error('Must specify dimension for matrix input');
    end
end
smArgs = varargin(desArgs);

% run SMOOTH over dimension requested
if tDim == 1
    nSteps = matSize(2);
else
    nSteps = matSize(1);
end
outMat = mat*NaN;

for iE = 1:nSteps
    if tDim == 1
        tIxNs = sub2ind(matSize, ...
                        1:matSize(1), ...
                        repmat(iE,[1 matSize(1)]));
    else
        tIxNs = sub2ind(matSize, ...
                        repmat(iE,[1 matSize(2)]), ...
                        1:matSize(2));
    end
    
    tV = smooth(mat(tIxNs), smArgs{:});
    outMat(tIxNs) = tV;
end

output = outMat;
