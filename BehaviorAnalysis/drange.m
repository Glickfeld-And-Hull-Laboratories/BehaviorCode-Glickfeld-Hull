function y = drange(x, dim)
%DRANGE (ps-utils): range of data along given dimension
%   Y = DRANGE(X, DIM)
%   DIM defaults to 1.  
% 
%  MH - http://github.com/histed/tools-mh

if nargin < 2, dim = 1; end
y = max(x,[],dim) - min(x,[],dim);
