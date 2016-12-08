function [outS] = behavDataExtractOneBlock(bs, block2N, side)
% given bs from getBehavDataForData, extract one set of params for block2
%
% histed 121004

outS = bs;

if bs.nBlock2Indices == 1 && bs.nside == 1
    outS = bs; 
    return
end


assert(block2N >= 1 && block2N <=2, 'b2N must be in range [1,2]');

matFields = { 'reactTimesByPower', 'reactTimeMean', 'reactTimeStd', 'reactTimeSEM', ...
    'percentsCorrect', 'nMiss', 'nCorr', 'nCorrPlusMiss', 'intensitiesC', 'percentsCorrectC', ...
    'nCorrPlusMissC'};
%cellExtractFields = { 'intensitiesC' };
cellExtractFields = {};

% do fields that are 2xN or Nx2
for iF = 1:length(matFields)
    tFN = matFields{iF};
    tF = bs.(tFN);
    desDim = size(tF);
    
    if length(desDim)<3
        bs.(tFN) = tF{block2N};
        tF = bs.(tFN);
        desDim = size(tF);
    end
    
    assert(~isempty(desDim), ['bug: field #' num2str(iF) ' not found']);
    
    if bs.nBlock2Indices > 1 && bs.nside == 1
        ind = find(desDim == 2);
        if ind == 1
        	outS.(tFN) = bs.(tFN)(block2N,:);
        else
            outS.(tFN) = bs.(tFN)(:,block2N);
        end
    elseif bs.nBlock2Indices == 1 && bs.nside > 1
        ind = find(desDim == 1);
        if ind == 1
            outS.(tFN) = bs.(tFN)(block2N,:,side);
        else
            outS.(tFN) = bs.(tFN)(:,block2N,side);
        end
    elseif bs.nBlock2Indices > 1 && bs.nside > 1 %check orientation of these matrices
        ind = find(desDim == 2);
        if ind == 1
            outS.(tFN) = bs.(tFN)(block2N,:,side);
        else
            outS.(tFN) = bs.(tFN)(:,block2N,side);
        end
    else
        error('bug: too many dims?');
    end
end

for iF = 1:length(cellExtractFields)
    tFN = cellExtractFields{iF};
    outS.(tFN) = bs.(tFN){block2N, side};       % extracts the cell contents

end


% special case a few
outS.block2Indices = bs.block2Indices(block2N);
outS.nBlock2Indices = 1;