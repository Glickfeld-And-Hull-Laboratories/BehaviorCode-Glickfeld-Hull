function y = nansum_dim(x,dim)
%NANSUM_DIM (ps-utils): sum value ignoring NaNs in any dimension
%   Y = NANSUM(X,DIM)
%   Returns the sum, ignoring NaNs. 
%   DIM specifies the dimension to work along.  Defaults to first
%   non-singleton dimension.
%
%   Notes: not needed for matlab 6 (R14)
%
%   MH - http://github.com/histed/tools-mh

x(isnan(x)) = 0;
y = sum(x,dim);

