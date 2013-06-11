function X = floorto(Xin,n)
%FLOORTO (ps-utils) take floor to n decimal places
%   FLOORTO(X,n) like ROUNDTO but takes the floor rather than rounding
%   e.g. FLOORTO(3.141592,4) returns 3.141500000..
%
%   See also CHOP, ROUNDTO, CEILTO
%
%  MH - http://github.com/histed/tools-mh

X = floor(Xin .* 10.^(n)) ./ 10.^(n);
