function outStr = sprintf_vector(pat, varargin)
%   outStr = SPRINTF_VECTOR(pat, inVect1, inVect2, ...)
%
%   Given a pattern and a vector, return a cellstr, with sprintf(pat,...) called on
%   each entry in the vector.
%
%   histed 111219
%
% MH - http://github.com/histed/tools-mh

nVects = length(varargin);
allLens = cellfun(@length, varargin);
vectLen = max(allLens);
assert(all(allLens == vectLen), 'all input vectors must have the same length');

vectC = cell(nVects, vectLen);
for iV = 1:nVects
    tV = varargin{iV};
    
    if isa(tV, 'numeric')
        tV = mat2cell_singleton(tV);
    elseif isa(tV, 'cell')
        % pass
    else
        error('not implemented');
    end
    if length(tV) == 1
        tV = repmat(tV, [vectLen, 1]);
    end
    
    vectC(iV,1:vectLen) = tV(:);
end

outStr = {};
% probably a more succinct way to do this w/ cellfun but this works
for iE = 1:vectLen
    outStr{iE} = sprintf(pat, vectC{:,iE});
end

