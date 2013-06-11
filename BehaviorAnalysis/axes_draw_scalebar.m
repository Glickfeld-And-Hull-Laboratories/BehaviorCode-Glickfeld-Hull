function sbH = axes_draw_scalebar(varargin)
%AXES_DRAW_SCALEBAR (mh-utils): convert tick-labeled axes into scalebar labeled
%   sbH = AXES_DRAW_SCALEBAR(varargin)
%
%   This function will turn off the axes and instead draw a scalebar with
%   the correct units
%
%   Defaults:
% defs = { 'XLength', chop(xTick(end)-xTick(end-1),1), ...
%          'YLength', chop(yTick(end)-yTick(end-1),1), ...
%          'CenterPosXY', [], ...
%          'Location', 'SouthEast', ...
%          'LineStyle', { 'LineWidth', 3, ...
%                         'Color', 'k' }, ...
%          'TextStyle', { }, ...
%          'XNumFmtStr', '%g', ...
%          'YNumFmtStr', '%g', ...         
%          'XLabel', '', ...
%          'YLabel', '', ...
%        };
%
%   specify bar location with 'Location'  ('NorthWest' etc.)
%       or 'CenterPosXY'  (in axes units, not normalized)
%
%   Example
%   figH = figure;
%   plot(rand(1,100));
%   lH = axes_draw_scalebar('XLabel', 'points', ...
%                           'YLabel', 'rand units', ...
%                           'Location', 'NorthWest');
%   
%
%  MH - http://github.com/histed/tools-mh

axH = gca;
figH = get(axH, 'Parent');
hold on;

xTick = get(axH, 'XTick');
yTick = get(axH, 'YTick');

defs = { 'XLength', chop(xTick(end)-xTick(end-1),1), ...
         'YLength', chop(yTick(end)-yTick(end-1),1), ...
         'CenterPosXY', [], ...%[xTick(end-1), yTick(1)], ...
         'Location', 'SouthEast', ...
         'LineStyle', { 'LineWidth', 3, ...
                        'Color', 'k' }, ...
         'TextStyle', { }, ...
         'XNumFmtStr', '%g', ...
         'YNumFmtStr', '%g', ...         
         'XLabel', '', ...
         'YLabel', '', ...
         'XLabelSide', 'bottom', ...
         'YLabelSide', 'left',
       };

uo = stropt2struct(stropt_defaults(defs, varargin));
oIn = stropt2struct(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set center position
if ~isempty(uo.CenterPosXY)
    assert(~isfield(oIn, 'Location'), ...
           'Cannot specify both CenterPosXY and Location');
    centerPosXY = uo.CenterPosXY;
else
    switch uo.Location
      case 'SouthWest'
        centerPosXY = [xTick(1) yTick(1)];
      case 'SouthEast'      
        centerPosXY = [xTick(end-1) yTick(1)];
      case 'NorthWest'
        centerPosXY = [xTick(1) yTick(end-1)];        
      case 'NorthEast'
        centerPosXY = [xTick(end-1) yTick(end-1)];    
      case 'South'
        centerPosXY = [(xTick(1)+xTick(end-1))/2 yTick(1)];
      case 'North'      
        centerPosXY = [(xTick(1)+xTick(end-1))/2 yTick(end-1)];
      case 'West'
        centerPosXY = [xTick(1) (yTick(1)+yTick(end-1))/2];        
      case 'East'
        centerPosXY = [xTick(end-1) (yTick(1)+yTick(end-1))/2];
        otherwise
            error('Invalid location value %s', uo.Location)
    end
end
    
% first turn off normal labels
set(axH, 'Visible', 'off');
% if set to manual, scale bar may be clipped 
% but we will live with this; otherwise manually set limits will be reverted
%         , ...
%          'XLimMode', 'auto', ...  
%          'YLimMode', 'auto');

set(figH, 'Color', 'w');

% draw line
xLen = uo.XLength;   % chop only if default value used
yLen = uo.YLength;
ptX = centerPosXY(1) + [0 0 xLen];
ptY = centerPosXY(2) + [yLen 0 0];
lH = plot(ptX, ptY);
if ~isempty(uo.LineStyle)
    set(lH, uo.LineStyle{:});
end


% label text
switch lower(uo.YLabelSide)
    case 'left'
        tYPos = centerPosXY + [-uo.XLength/10 uo.YLength/4 ];
        yHAlign = 'right';
    case 'right'
        tYPos = centerPosXY + [+uo.XLength/7 uo.YLength/2 ];
        yHAlign = 'left';
    otherwise
        error('Invalid value for YLabelSide: %s', uo.YLabelSide);
end
switch lower(uo.XLabelSide)
    case 'bottom'
        tXPos = centerPosXY + [uo.XLength/2  -uo.YLength/10];
        xVAlign = 'top';
    case 'top'
        tXPos = centerPosXY + [uo.XLength/2  +uo.YLength/7];
        xVAlign = 'bottom';
    otherwise
        error('Invalid value for YLabelSide: %s', uo.YLabelSide);
end





xStr = [sprintf(uo.XNumFmtStr,xLen) ' ' uo.XLabel];
yStr = [sprintf(uo.YNumFmtStr,yLen) ' ' uo.YLabel];

if uo.XLength > 0
    tH(1) = text(tXPos(1), tXPos(2), xStr, ...
                 'VerticalAlignment', xVAlign, ...
                 'HorizontalAlignment', 'center');
end

if uo.YLength > 0
    tH(2) = text(tYPos(1), tYPos(2), yStr, ...
                 'VerticalAlignment', 'bottom', ...
                 'HorizontalAlignment', yHAlign);
end
    
if ~isempty(uo.TextStyle)
    set(tH, uo.TextStyle{:});
end


% export handles
sbH = cat(2,lH, tH);
set(sbH, 'Tag', 'axes_draw_scalebar');
