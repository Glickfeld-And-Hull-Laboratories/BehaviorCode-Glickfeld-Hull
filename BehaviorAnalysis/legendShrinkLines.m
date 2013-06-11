function legendShrinkLines(legH, pctShrink)
%% function legendShrinkLines(legH, pctShrink)
% must come after all calls to change lines in figure, else the new length
% will be overwritten
% - pctShrink == 100 means do not change
% 
% histed 121008

%% test
lH = findobj(legH, 'Type', 'line');

xdC = get(lH, 'XData');  % all lines used for drawing have len 2
legLineIx = cellfun(@(x) length(x)==2, xdC);
lH = lH(legLineIx);
nL = length(lH);

for iL=1:nL
    xd =get(lH(iL), 'XData');
    tLen = range(xd);
    newLen = tLen*pctShrink/100;
    newXd = [xd(2)-newLen xd(2)];
    set(lH(iL), 'XData', newXd);
end

