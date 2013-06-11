function y = rectify(x)
%RECTIFY (ps-utils): if input < 0, make it equal 0, if positive, return value
%  Y = RECTIFY(X)
% 
%  MH - http://github.com/histed/tools-mh

x(x<0) = 0;
y = x;