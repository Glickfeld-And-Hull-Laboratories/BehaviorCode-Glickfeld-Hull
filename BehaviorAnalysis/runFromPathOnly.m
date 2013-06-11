function varargout=runFromPathOnly(fcnName, varargin)
%runFromPathOnly: execute function in path, skipping any in current directory
% 
%   [fcnArgsOut] = runFromPathNotCurrDir(fcnName, fcnArgsIn)
%       fcnName is a string, not a function handle
%
%       fcn is not evaluated in the caller's context; the only variables it can use is those 
%           passed as arguments.  This is equivalent to using feval(fcnName, args).
%
% histed 121225

assert(ischar(fcnName), 'fcnName must be a string');

p = which(fcnName, '-all');
currDir = pwd;
[currDirFromPath,~,~] = fileparts(p{1});
if ~strcmp(currDir, currDirFromPath)
    warning('%s.m not found in current dir, running from path', fcnName);
    [defaultDir,~,~] = fileparts(p{1});  % use first element of path
else
    [defaultDir,~,~] = fileparts(p{2});  % use first element of path that's not curr dir
end

try
    cd(defaultDir);
    nOut = nargout;
    if nOut == 0, nOut = 1; end
    [outC{1:nOut}] = feval(fcnName, varargin{:});
    cd(currDir);
    varargout(1:nOut) = outC(1:nOut);
    return
catch me
    % cleanup in case of errors
    cd(currDir);
    rethrow me;
end
