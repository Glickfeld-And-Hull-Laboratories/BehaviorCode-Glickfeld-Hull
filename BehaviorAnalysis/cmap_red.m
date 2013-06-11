function scaledMap = cmap_red(m)
%CMAP_RED (ps-utils): colormap linear on green values
%
%  MH - http://github.com/histed/tools-mh

if nargin < 1, m = size(get(gcf, 'colormap'), 1); end

scaledMap = [0:(m-1); zeros(1,m); zeros(1,m)]'./(m-1);

