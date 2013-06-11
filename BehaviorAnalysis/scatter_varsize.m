function outH = scatter_varsize(x, y, scalesize, extraparms)
%SCATTER_VARSIZE (ps-utils): scatter plot with replications 
%   SCATTER_VARSIZE(X, Y, SCALESIZE, EXTRAPARMS)
%
%   See also: SCATTER
%
%  MH - http://github.com/histed/tools-mh

if nargin < 3 || isempty(scalesize), scalesize = 3; end
if nargin < 4, extraparms = {}; end
    
x=colvect(x);
y=colvect(y);

[unq unqix] = unique([x, y], 'rows');
sizes = ones(length(unqix),1);

% probably can be optimized by using histc
for i=1:length(unqix)
    sizes(i) = sum( (x == unq(i,1)) & (y == unq(i,2)) );
end

pointH = scatter(unq(:,1), unq(:,2), sizes .* scalesize, extraparms{:});
if nargout > 0, outH = pointH; end
    
