function outS = concatenateDataBlocks(blockC)
%% combining multiple mat files

    fNames = fieldnames(blockC(1));
    nF = length(fNames);
    outS = struct([]);
    trs = blockC(1).trialSinceReset;
    for iF = 1:nF
        fN = cell2mat(fNames(iF,:));
        outS(1).(fN) = blockC(1).(fN);
        if size(blockC(1).(fN),2) == trs;
            for iblock = 2:size(blockC,2)
                outS.(fN) = [outS.(fN) blockC(iblock).(fN)];
            end
        else
            outS.(fN) = blockC(1).(fN);
        end
    end