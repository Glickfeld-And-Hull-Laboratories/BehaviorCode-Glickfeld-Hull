function varargout = call_function_in_dir(fName, workingDir, varargin)
%CALL_FUNCTION_IN_DIR (ps-utils): call a function in a given working dir
%   varargout = CALL_FUNCTION_IN_DIR(fHandle, workingDir, varargin)
%
%   Handles cleanup 
% 
%  MH - http://github.com/histed/tools-mh

oldDir = pwd;

% change to new dir
cd(workingDir);

fHandle = str2func(fName);
nFuncArgOut = nargout(fHandle);

try
    [varargout{1:nFuncArgOut}] = feval(fHandle, varargin{:});
catch
    % error in execution...

    % restore directory
    cd(oldDir);    
    rethrow(lasterror);
end


% success
cd(oldDir);




