function outPath = fullfileMH(varargin)
%fullfileMH (ps-utils): concatenate paths, using '/' always as filesep
%
%   Windows MATLAB understands '/' as a filesep, and sprintf gets confused
%   by '\', so use '/'
%
%  MH - http://github.com/histed/tools-mh

% first strip trailing slashes
varargin = regexprep(varargin, '(.*)/$', '$1');

% concat
nPathsIn = length(varargin);
newOut = cell(1,nPathsIn*2-1);
newOut(1:2:end) = varargin(:);
[newOut{2:2:end}] = deal('/');
outPath = cat(2,newOut{:});
