function outS = concatenateStructuresLG(blockC)
%% combining multiple mat files

    fNames = fieldnames(blockC(1));
    nF = length(fNames);
    outS = struct;
    for iF = 1:nF
        fN = cell2mat(fNames(iF,:));
        outS(1).(fN) = blockC(1).(fN);
        for iblock = 2:size(blockC,2)
            outS.(fN) = [outS.(fN) blockC(iblock).(fN)];
        end
    end