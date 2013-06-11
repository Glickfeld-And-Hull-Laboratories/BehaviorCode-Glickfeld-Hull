function edit (fName)
%EDIT (tools-mh): open a matlab file in your favorite editor
%
% MH 100826 - created

persistent isOpened

if isempty(isOpened)
    isOpened = false;
end

fullName = which(fName);
if isempty(fullName)
    error('''%s'' not found on path', fName);
elseif strcmp(fullName(1:8), 'built-in')
    error('''%s'' is a Matlab built-in', fName);
end


% $$$ if isOpened
% $$$     tEditor = '/Applications/Aquamacs.app/Contents/MacOS/bin/emacsclient';
% $$$ else
% $$$     tEditor = '/Applications/Aquamacs.app/Contents/MacOS/Aquamacs';
% $$$ end

if ~isOpened
    unix('/Applications/Aquamacs.app/Contents/MacOS/Aquamacs -f server-start &');
end

tEditor = '/Applications/Aquamacs.app/Contents/MacOS/bin/emacsclient';
unix(sprintf('%s %s &', tEditor, fullName));
isOpened = true;
