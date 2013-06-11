function handles = vert_shade_gauss (mu, sigma, maxintensity, yminmax)
%VERT_SHADE (ps-utils): draw vertical gaussian-shaded square gradient 
%   HANDLES = VERT_SHADE_GAUSS (MU, SIGMA, MAXINTENSITY, YMINMAX)
%   
%   See also VERT_SHADE.
%
%  MH - http://github.com/histed/tools-mh

if nargin < 4, yminmax = []; end

x0 = mu - 3*sigma;
x1 = mu;
x2 = mu + 3*sigma;

increm = (x2-x0) / 1000;

xp = [x0:increm:x1, x1:increm:x2];

% construct gaussian color vector
nx = length(xp);
gaussc = normpdf(xp(1:end-1), mu, sigma);
gaussc = gaussc / max(gaussc) * maxintensity  ;
hsvc = ones(length(xp)-1, 3);
hsvc(:,2) = gaussc';
colors = ones(1,length(xp)-1, 3);
colors(1,:,:) = hsv2rgb(hsvc);

vert_shade(xp, colors, yminmax);

