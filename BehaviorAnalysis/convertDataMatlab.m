function dsOut = convertDataMatlab(ds)
% Sanitize the matlab data structure: convert all fields to structs, etc.

% MH 130118

%% process vectors in input
%vFields = { 'tLaserPowerMw', 'tGratingContrast', 'tTotalReqHoldTimeMs', ...
%    'tTotalReqHoldTimeMs', 'holdTimesMs', 'reactTimesMs' };
fNames = fieldnames(ds);
for iV=1:length(fNames)
    tFN = fNames{iV};
    if iscell(ds.(tFN))
        try
            dsOut.(tFN) = celleqel2mat_padded(ds.(tFN), NaN, 'double');
        catch
            dsOut.(tFN) = ds.(tFN);
        end
    else
        dsOut.(tFN) = ds.(tFN);
    end
end

