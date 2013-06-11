function cellout = ensure_cellstr (string)
%ensure_cellstr (ps-utils): make sure a string is a cell matrix of strings
%
% USAGE
%  function cellout = ensure_cellstr (string)
% 
% If the input is a cell matrix it's left alone.
% If it's a regular char matrix it's made into a cell matrix of strings.
%
%   MH - http://github.com/histed/tools-mh

if iscellstr (string)
        cellout = string;
        return
elseif ischar (string)
        cellout = cellstr (string);
        cellout = deblank (cellout);
        return
else
        error ('Not a string of either sort!');
end


