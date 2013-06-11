function trialOutcomeCell = subComputeTrialOutcomeCell(ds)
%  trialOutcomeCell = subComputeTrialOutcomeCell(ds)
%  ds is output of convertDataMatlab/Serial
% 
% MH 130117

rc = behavConstsHADC8;

foundOutcomes = intersect(rc.trialOutcomes, fieldnames(ds));

nO = length(foundOutcomes);
if nO == 0
    error('No outcome vectors found in data struct');
end

for iO = 1:nO
    tOutcomeStr = foundOutcomes{iO};
    outcomeIxC{iO} = diff([0 ds.(tOutcomeStr)]);
end
outcomeMat = cat(1,outcomeIxC{:});

[nOutcomes nTrials] = size(outcomeMat);
totOut = sum(outcomeMat,1);
assert(all(totOut>0), 'zero outcomes per trial?!');

    

trialOutcomeCell={};
for iT = 1:nTrials
    outcomeN = find(outcomeMat(:,iT)>=1);
    if length(outcomeN) > 1
        warning('multiple outcomes per trial!');
    end
    trialOutcomeCell{iT} = foundOutcomes{outcomeN};
end

    
    
