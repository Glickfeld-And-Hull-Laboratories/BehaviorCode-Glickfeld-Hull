function [axH, figH] = figs_of_subplots(spSize1, spSize2, spN);
%
%   for iC = 1:maxN
%       axH = figs_of_subplots(3,3,iC);
%   end
%
% histed 111201

totSp = spSize1*spSize2;

nextSpN = mod(spN-1, totSp)+1;

if nextSpN == 1;
    fig1H = figure;
end

axH = subplot(spSize1, spSize2, nextSpN);

if nargout > 1
    figH = fig1H;
end
