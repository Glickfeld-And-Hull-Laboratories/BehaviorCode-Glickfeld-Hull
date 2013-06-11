function titlestr = make_title (extraText, mName)
%MAKE_TITLE (stacks_mh): return a string to be used as a plot title
%
%  (note: simplified version of MAKE_PLOT_TITLE (posit)
%
%   MH - http://github.com/histed/tools-mh

if nargin < 2 || isempty(mName)
    mName = caller_mfilename(1, true);
end
mName = strrep(mName, '_', '\_');

titlestr = sprintf('\\bf%s\\rm: %s', mName, extraText);

