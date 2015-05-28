function clipclip(figH, tSize, fFormat, boundingBox)
%CLIPCLIP (tools): copy current figure to a snapshot on disk for ppt etc.
%   clipclip(figH, tSize, fFormat)
% 
%   export path is hardcoded below
%
%   notes: figures with transparency don't work with this, printing
%   doesn't work (matlab bug).  Use clipclippix to grab a screenshot.
%
%   See also CLIPCLIPPIX
%
%  MH - http://github.com/histed/tools-mh

if nargin < 1 || isempty(figH), figH = gcf; end
if nargin < 2 || isempty(tSize), tSize = 6*[1 0.75]; end
if nargin < 3 || isempty(fFormat), fFormat = 'pdf'; end


%% set paths
% file export path
%fPath = '~/shared/snapshots';
%fPath = 'i:/users/histed/snapshots';
%dirs = directories;
%assert(isfield(dirs, 'toolsSnapshots'), ....
%       'needed dir entry is missing, edit directories.m');

nFigs = length(figH);
allNameC = {};
for iFig = 1:nFigs
    tFig = figH(iFig);

    fName = [mfilename '-' datestr(now, 'yymmdd-HHMM_SS_FFF')];
    allNameC{iFig} = fName;
    
    fullF = fullfileMH('C:\Users\andrew\Desktop', fName);
    exportfig_print(tFig, fullF, ...
                    'FileFormat', fFormat, ...
                    'Size', tSize, ...
                    'Renderer', 'painters');
    
end
if length(allNameC) > 1
    if any(strcmp(allNameC{1}, allNameC(2:end)))
        error('clipclip working too fast, need to change filename to avoid collison');
    end
end    





%unix(sprintf('convert %s %s', [fullF '.png'], [fullF '.wmf']));
