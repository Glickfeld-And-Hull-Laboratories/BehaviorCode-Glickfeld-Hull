function structArray = structincell2structarray(cellArrayIn, fillVal)
%STRUCTINCELL2STRUCTARRAY (ps-utils): concat structures inside cell array
%   structArray = STRUCTINCELL2STRUCTARRAY(cellArrayIn, fillVal)
% 
%   cellArrayIn must have a single structure in each cell.  
%   If a cell contains NaN (the
%   singleton numeric matrix), all fields of the struct for that cell are set
%   to fillVal.
%   fillVal defaults to [], NaN is also useful
%
%   If fields are different, corresponding entries have empty field values
%
%  MH - http://github.com/histed/tools-mh

fillVal = [];

% find 1x1 NaN numeric matrices in CellArrayIn 
nonStructIx = ~cellfun('isclass', cellArrayIn, 'struct');
nonStructVect = cat(1,cellArrayIn{find(nonStructIx)});
assert(all(isnumeric(nonStructVect)) && all(isnan(nonStructVect)), ...
       'Expect all non-struct values in cell array to be NaNs !');
nanIx = nonStructIx;


% replace singleton nans with the filler struct
%filledCell = cellArrayIn;
%filledCell(nanIx) = { fillerStruct };

% find field names for each
structCell = cellArrayIn(~nanIx);

% find added field names
fN1 = fieldnames(structCell{1}); % optimization
fNameAdd = cellfun(@(x) setdiff(fieldnames(x),fN1), ...
                 structCell, 'UniformOutput', false);
% subtracted
fNameSub = cellfun(@(x) setdiff(fN1, fieldnames(x)), ...
                   structCell, 'UniformOutput', false);
addNames = unique(cat(1,fNameAdd{:}));
subNames = unique(cat(1,fNameSub{:}));
for iC=1:length(structCell)
    tExtras = setdiff(addNames, fNameAdd{iC});
    tMissing = fNameSub{iC};
    tS = structCell{iC};
    tList = union(tExtras,tMissing);
    for iE = 1:length(tList)
        tS.(tList{iE}) = fillVal; 
    end
    structCell{iC} = orderfields(tS);
end


sortC = cellfun(@(x) x, structCell);  % dummy fn to return the structure

% make a filler struct where all fields are NaN
structFNames = fieldnames(sortC);
nFNames = length(structFNames);
for iF=1:nFNames
    tFName = structFNames{iF};
    fillerStruct.(tFName) = fillVal;
end

% replace blanks with fill struct
outS(1) = sortC(1);
outS(find(~nanIx)) = sortC(:);
outS(nanIx) = deal(fillerStruct);
structArray = outS;


