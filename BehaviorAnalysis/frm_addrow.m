function xdOut = frm_addrow(xd, afterRowN)
%FRM_ADDROW: (tools-mh): add a blank row to a frame structure
%
%   xdOut = frm_addrow(xd, afterRowN)
%
% histed 120626

if nargin < 2, afterRowN = []; end

if ~isempty(afterRowN), error('afterRowN not implemented now - fix this'); end

allFields = setdiff(fieldnames(xd), { 'colNames', 'nCols', 'nRows' });
nFields = length(allFields);
nextRowN = xd.nRows+1;

xdOut = xd;
for iF = 1:nFields
    tFN = allFields{iF};
    tFV = xd.(tFN);
    if isa(tFV, 'double')
        xdOut.(tFN)(nextRowN) = NaN;
    elseif iscell(tFV)
        xdOut.(tFN){nextRowN} = [];
    else
        error('bug: unknown field type');
    end
end
xdOut.nRows = nextRowN;

