function axH = subplotFromXY(maxXYVals, xyVals)
%
%   Call this as axH = subplotFromXY([maxX maxY], [x,y])   
%
% histed 111014

tX = xyVals(1);
tY = xyVals(2);
if any(maxXYVals <= 1)
    error('invalid xyVals');
end

maxX = maxXYVals(1);
maxY = maxXYVals(2);
if any(maxXYVals <= 1)
    error('invalid maxXYVals');
end

% very simple code but tricky to get right (i.e. I always screw this up), so
% wrap it in a function
p = sub2ind([maxX, maxY], tX, tY);
axH = subplot(maxY, maxX, p);


% subplot(m,n,p) - m is rows/y; n is cols/x, but p runs along the cols
