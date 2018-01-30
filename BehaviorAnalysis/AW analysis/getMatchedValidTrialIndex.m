function validIndSamples = getMatchedValidTrialIndex(...
    tValidTargets,tInvalidTargets)

tValidTargets = round(double(tValidTargets), 2, 'significant');
valTargets = unique(tValidTargets);
tInvTargets = round(double(tInvalidTargets), 2, 'significant');
invTargets = unique(tInvTargets);
nInvTargetType = hist(tInvTargets,invTargets);
invTargetInd = ismember(valTargets,invTargets);

valTargets = valTargets(2:end);
invTargets = invTargets(2:end);
nInv = length(invTargets);
invTargetInd = invTargetInd(2:end);
nInvTargetType = nInvTargetType(2:end);

nBoot = 10;
if nInv == 1
    valIndSample = {find(tValidTargets == valTargets(invTargetInd))};
else
    valIndSample = cell(1,nInv);
    for i = 1:nInv
        if ~any(valTargets == invTargets(i))
            continue
        else
            ind = find(tValidTargets == valTargets(valTargets == invTargets(i)));
            if length(ind) < 2
                ind = [ind ind]
            end
            visIndBoot = cell(nBoot,1);
            for iboot = 1:nBoot
                visIndBoot{iboot,1} = randsample(ind,nInvTargetType(i));
            end
            valIndSample{i} = cell2mat(visIndBoot');
        end
    end
end

validIndSamples = valIndSample;
end

