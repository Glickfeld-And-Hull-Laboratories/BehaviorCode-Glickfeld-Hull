function [matchValTrials nMatch_val] = matchTrials(invVisTargets,invTars,invTrials,visTargets);

matchValTars = ismember(visTargets,invTars);
nTar_inv = histcounts(invVisTargets(invTrials),length(invTars));
nTar_val = histcounts(visTargets(matchValTars),length(invTars));
if length(invTars) == 2;
    invRatio = nTar_inv(1)/nTar_inv(2);
    valRatio = nTar_val(1)/nTar_val(2);
    nMatch_val = zeros(1,2);
    
    if invRatio > valRatio
        nMatch_val(1) = nTar_val(1);
        nMatch_val(2) = round(nTar_val(1) * (1/invRatio));
    else
        nMatch_val(2) = nTar_val(2);
        nMatch_val(1) = round(nTar_val(2) * invRatio);
    end
        
    matchValTrials_type = cell(1,2);
    for i = 1:2
        t = visTargets == invTars(i);
        ind = randsample(sum(t),nMatch_val(i));
        tr = find(t);
        matchValTrials_type{i} = tr(ind);
    end
    matchValTrials = sort(cell2mat(matchValTrials_type));
elseif length(invTars) == 1
    matchValTrials = find(matchValTars);
else
    error('length > 2')    
end
end
