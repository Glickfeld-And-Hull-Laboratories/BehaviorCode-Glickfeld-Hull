function X = ceilto(Xin,n)
%CEILTO (ps-utils) take ceil to n decimal places
%   CEILTO(X,n) like ROUNDTO but takes the ceil rather than rounding
%   e.g. CEILTO(3.141592,4) returns 3.141500000..
%
%   See also CHOP, ROUNDTO, FLOORTO
%
%  MH - http://github.com/histed/tools-mh

X = ceil(Xin .* 10.^(n)) ./ 10.^(n);
