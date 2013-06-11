function [ebH whiskerH] = basicerrorbar(lH, upperLength, lowerLength, whiskerLen)
%BASICERRORBAR (mh): construct errorbars with only line objects
%    [ebH whiskerH] = basicerrorbar(lineH, upperLength, lowerLength, whiskerLen)
%  
%    Main advantage of this over errorbar is that only line objects are used,
%    so it can be easily modified
%
%    whiskerLen defaults to 0, meaning do not draw whiskers.  Length is in X
%    axis units.
%
% histed 120529 brought in from plotOnlineHist.m

if nargin < 3 || isempty(lowerLength), lowerLength = upperLength; end
if nargin < 4, whiskerLen = 0; end

nB = length(lH);
lbDims = size(upperLength);
desND = find(lbDims ~= nB);
if desND == 1
    upperLength = upperLength';
    lowerLength = lowerLength';
end


isLog = strcmp(get(gca, 'XScale'), 'log');

for iB=1:nB
    xv = get(lH(iB), 'XData');
    yv = get(lH(iB), 'YData');
    y0 = yv-lowerLength(iB,:);
    y1 = yv+upperLength(iB,:);

    xvals = cat(1, xv, xv, xv);
    yvals = cat(1, y0, y1, yv*NaN);

    hold on;
    ebH(iB,1) = plot(xvals(:)', yvals(:)');

    if whiskerLen > 0
        wl = whiskerLen/2;
        if isLog
            wl = wl+1;
            wxvals = cat(1, xv./wl, xv.*wl, xv, xv./wl, xv.*wl, xv);
        else % linear
            wxvals = cat(1, xv-wl, xv+wl, xv, xv-wl, xv+wl, xv);
        end
        wyvals = cat(1, y0, y0, y0*NaN, y1, y1, y1*NaN);
        whiskerH(iB,1) = plot(wxvals(:)', wyvals(:)');
    end
end
propName = 'Color';
if length(lH) > 1, propName = {propName}; end
set(ebH, propName, get(lH, 'Color'));
set(ebH, 'Tag', 'basicerrorbar_stem');

if whiskerLen > 0
    set(whiskerH, propName, get(lH, 'Color'));
    set(whiskerH, 'Tag', 'basicerrorbar_whisker');
end

