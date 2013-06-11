function pH = lim_lines(posVect, alongAxesStr, varargin)
%
%   pH = lim_lines(posVect, alongAxesStr, lineprops)
%
%   Like vert_lines but works on x or y axes
%   For vertical lines, alongAxesStr = 'Y'
%
%   for update:
%   lim_lines('update', figH)  - call on the overlying figure
%
% histed 111016

%% do update if requested
if nargin > 1 && ischar(posVect) && strcmpi(posVect, 'update')
    figH = alongAxesStr;
    if ~strcmp(get(figH, 'Type'), 'figure')
        error('pass in a figure handle as 2nd parameter on update');
    end
    subDoUpdate(figH, 'X');
    subDoUpdate(figH, 'Y');

    return
end

%% process normally if not update
if nargin < 2 || isempty(alongAxesStr)
    alongAxesStr = 'X'; ...
end
if nargin < 3, 
    lineProps = {'k'}; 
else
    lineProps = varargin;
end

%% check errors
if ischar(alongAxesStr) & ~any(strcmpi({'x', 'y'}, alongAxesStr))
    error('invalid value for alongAxesStr: %s', alongAxesStr);
end




%%
tagStr = sprintf('%sLim', upper(alongAxesStr));

tLim = get(gca, tagStr);
switch upper(alongAxesStr)
  case 'X'
    pH = plot((ones(size(posVect(:)))*tLim)', ...
              (posVect(:)*[1 1])', ...
              lineProps{:});
  case 'Y'
    pH = plot((posVect(:)*[1 1])', ...
              (ones(size(posVect(:)))*tLim)', ...
              lineProps{:});
end


set(pH, 'Tag', sprintf('%sLimLines', alongAxesStr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subDoUpdate(figH, alongAxesStr);    

lList = findobj(gcf, 'Tag', sprintf('%sLimLines', alongAxesStr));
nL = size(lList,1);
for iL = 1:nL
    tH = lList(iL);
    axH = get(tH, 'Parent');
    tLim = get(axH, [alongAxesStr 'Lim']);
    set(tH, [alongAxesStr 'Data'], tLim);
end

