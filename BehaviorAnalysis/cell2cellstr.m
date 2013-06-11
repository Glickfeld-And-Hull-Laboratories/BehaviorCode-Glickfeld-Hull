function cellstrOut = cell2cellstr(cellIn)
%CELL2CELLSTR (ps-utils): recursively convert all cells of a cell arr to str
%   CELLSTROUT = CELL2CELLSTR(CELLIN)
%   
%  MH - http://github.com/histed/tools-mh

assert(iscell(cellIn), 'Input must be a cell array');

cellstrOut = cell(size(cellIn));

nCells = prod(size(cellIn));
for i=1:nCells
    tEl = cellIn{i};
    if isnumeric(tEl)
        cellstrOut{i} = mat2str(tEl);
    elseif iscell(tEl)
        cellstrOut{i} = cell2cellstr(tEl);
    else
        error('Elements of cellIn must be numeric or cell arrays');
    end
end

    
