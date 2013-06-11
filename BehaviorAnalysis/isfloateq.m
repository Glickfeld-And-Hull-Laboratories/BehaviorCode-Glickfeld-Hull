function bool = isfloateq(float1, float2, epsilon)
%ISFLOATEQ (ps-utils): are two floating point numbers equal?
%   BOOL = ISFLOATEQ(FLOAT1, FLOAT2)
%   or BOOL = ISFLOATEQ(FLOAT1, FLOAT2, EPSILON)
%   returns 1 if arguments differ by less than epsilon
%
%   EPSILON defaults to EPS*1e2, about 2.2e-14 for 8-byte doubles
%   see help EPS).
%
%   The goal is to for this to be interchangable with isequal for floats
%
%  MH - http://github.com/histed/tools-mh
if nargin < 3, epsilon = eps*1e2; end

bool = ( abs(float1 - float2) < epsilon );

