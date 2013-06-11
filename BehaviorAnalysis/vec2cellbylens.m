function vcell = vec2cellbylens(invec, vecLens)
%VEC2CELL (tools-mh): convert ragged vector into a cell vector
%   VCELL = VEC2CELLBYLENS(invec, vecLens)
%   This function takes a ragged vector stored as a vector and a set of lengths
%   and breaks the elements up into entries of a cell vector
%
%   See also: VEC2PADDED, VEC2PADBYCOL
%
%  MH - http://github.com/histed/tools-mh

nEls = length(vecLens);
vcell = cell([1 nEls]);


breakNs = [1 cumsum(vecLens)+1];
endNs = breakNs(2:end)-1;
startNs = breakNs(1:end-1);

for iL = 1:nEls
    tL = vecLens(iL);
    vcell{iL} = invec(startNs(iL):endNs(iL));
end

    
    
