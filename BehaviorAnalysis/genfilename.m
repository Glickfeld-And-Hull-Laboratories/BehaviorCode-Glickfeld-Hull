function outName = genfilename(inName)
%GENFILENAME (mh-tools): sanitize a string to be used as a filename
%
%   String must not have an extension.
%
% histed 111004

% $$$ [inPath inName inExt] = fileparts(inName);
% $$$ 
% $$$ if ~isempty(inPath)
% $$$     error('cannot specify a path');
% $$$     % could sanitize paths too but we do not now
% $$$ end
% the above fails because it uses periods to find extensions

v1 = regexprep(inName, '[^A-Z^a-z^0-9^\-^_]', '_');
outName = v1;


    
