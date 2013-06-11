function outFigH = dupfig(figH)
%DUPFIG (histed): make a new figure that is a copy of the current figure
%
%   outFigH = DUPFIG(figH)
%
% 110930 histed

if nargin < 1, figH = gcf; end

fig2H = figure;
copyobj(get(figH, 'Children'), fig2H);

if nargout > 0, outFigH = fig2H; end

