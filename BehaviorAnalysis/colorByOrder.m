function [tColor,nearestCmapNum] = colorByOrder(colorNum, axesH, colorOrder)
%COLORBYORDER (ps-utils): get next color from axes ColorOrder
%   [tColor,nearestCmapNum] = colorByOrder(colorNum, axesH, colormap)
%
%   Handles cycling through colors properly.
%
%   Can specify your own colormap/colorOrder as third argument
%
%  MH - http://github.com/histed/tools-mh

if nargin < 2, axesH = gca; end
if nargin < 3, colorOrder = []; end

if isempty(colorOrder)
    colorOrder = get(axesH, 'ColorOrder');
end

nColors = size(colorOrder, 1);

colorIx = mod(colorNum, nColors);
if colorIx == 0; colorIx = nColors; end

tColor = colorOrder(colorIx, :);

% get nearest cmap num if requested
figH = get(gca, 'Parent');
cmap = get(figH, 'Colormap');
nColorsInMap = size(cmap,1);

cmapDists = sum((cmap - repmat(tColor,nColorsInMap,1)) .^ 2, ...
                2);
[crap nearestCmapNum] = min(cmapDists);


    
