function outLim = equalize_lims(axH, alongAxesStr)
%EQUALIZE_LIMS (tools-mh): sets axes to same limits: [min max] of all
%
%   outLim = equalize_lims(axH, alongAxesStr)
%
%   Simple function but it saves me a little time debugging cell output of
%   get() etc. each time I use it
%
% histed 120220

%% args
if nargin < 2 || isempty(alongAxesStr)
    error('Must specific alongAxesStr');
end

%% check errors
if ~ischar(alongAxesStr) || ~any(strcmpi({'x', 'y'}, alongAxesStr))
    error('invalid value for alongAxesStr: %s', alongAxesStr);
end

%% do it
tagStr = sprintf('%sLim', upper(alongAxesStr));

assert(length(axH) > 1, 'error: only useful on multiple axes');
allLimC = get(axH, tagStr);
allLim = cat(1, allLimC{:});
newLim = [min(allLim(:,1)) max(allLim(:,2))];
set(axH, tagStr, newLim);

if nargout > 0
    outLim = newLim;
end
