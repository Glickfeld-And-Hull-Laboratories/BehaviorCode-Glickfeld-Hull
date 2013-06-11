function plot_distrib (inmat, flags)
%PLOT_DISTRIB (ps-utils): plot the mean and std dev of a matrix
%   PLOT_DISTRIB (INMAT, FLAGS)
%   inmat is a matrix whose columns are the vectors of interest.  
%   This function calculates the means and some measure of variance 
%   columnwise
%
%   flags:
%     'std': use std dev as the variance measure
%
%  MH - http://github.com/histed/tools-mh

if nargin < 2, flags = {'std'}; end

m = mean(inmat,2);

[stddev] = process_cell_flags({'std'}, flags);

if stddev
        vmlo = std(inmat, 0, 2);
        vmhi = vmlo;
end

plot (m);
hold on
plot (m-vmlo, ':');
plot (m+vmhi, ':');

