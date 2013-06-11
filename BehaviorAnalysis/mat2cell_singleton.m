function cellout = mat2cell_singleton (mat)
% mat2cell_singleton (ps-utils): mat2cell but cell has singleton entries
%
%  111219: Deprecated: use NUM2CELL.  This function was first created for a
%  very old version of matlab.
%
%  function cellout = mat2cell_singleton (mat)
%
% Same as:
%   mat2cell(mat, 1, ones (size(mat,1), size(mat,2)));
%   or
%   num2cell(mat)
%
%   MH - http://github.com/histed/tools-mh

%cellout = mat2cell(mat, ones(size(mat,1),1), ones(size(mat,2),1));
cellout = num2cell(mat);
