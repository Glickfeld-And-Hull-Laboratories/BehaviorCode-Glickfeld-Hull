function outCell = cellstr_cat(varargin)
%CELLSTR_CAT (ps-utils): take cellstr and char strings, return a single cellstr
%  OUTCELL = CELLSTR_CAT(VARARGIN)
%
%  MH - http://github.com/histed/tools-mh

% a lot of hairy work here to avoid looping.

tempCell = varargin;

charIx = cellfun('isclass', varargin, 'char');
nChars = sum(charIx);
allChars = varargin(charIx);
% wrap each string inside its own cell array
if ~isempty(allChars)
    buriedChars = mat2cell(allChars, 1, ones(1,nChars));
    % stick each wrapped char array back into tempCell
    tempCell(charIx) = buriedChars(:);
end
outCell=cat(2,tempCell{:});



