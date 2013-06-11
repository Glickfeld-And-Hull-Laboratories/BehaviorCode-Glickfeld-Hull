function scaledMap = cmap_blue(m)
%CMAP_GREEN (ps-utils): colormap linear on green values
%
%  MH - http://github.com/histed/tools-mh

if nargin < 1, m = size(get(gcf, 'colormap'), 1); end

scaledMap = [zeros(1,m); zeros(1,m); 0:(m-1)]'./(m-1);

