function [outA, outB, fnDiffA, fnDiffB] = fieldComp(dataA, dataB);

% finds differences in fields within two structures
% outA reports 1 for all fields in dataA that are also in dataB
% outB reports 1 for all fields in dataB that are also in dataA
% fnDiff is a cell array of fields that differ between structures

fNamesA = fieldnames(dataA);
nFA = length(fNamesA);

fNamesB = fieldnames(dataB);
nFB = length(fNamesB);

outA = zeros(1,nFA);
outB = zeros(1,nFB);
fnDiffA = [];
start = 1;
for iF = 1:nFA
    outA(:,iF) = isfield(dataB, fNamesA(iF));
    if outA(:,iF)== 0
        fnDiffA{start} = fNamesA(iF);
        start = start+1;
    end
end

fnDiffB = [];
start = 1;
for iF = 1:nFB
    outB(:,iF) = isfield(dataA, fNamesB(iF));
    if outB(:,iF)== 0
        fnDiffB{start} = fNamesB(iF);
        start = start+1;
    end
end