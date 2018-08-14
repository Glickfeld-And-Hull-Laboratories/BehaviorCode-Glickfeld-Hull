function outS = concatenateDataBlocks(blockC)
%% combining multiple mat files

    fNames = fieldnames(blockC(1));
    nF = length(fNames);
    outS = struct;
    
    for iblock = 1:size(blockC,2)
        if isfield(blockC, 'trialOutcomeCell')
            trs = length(blockC(iblock).trialOutcomeCell);
        elseif isfield(blockC, 'tGratingContrast')
            trs = length(blockC(iblock).tGratingContrast);
        elseif isfield(blockC, 'tStimOneGratingContrast')
            trs = length(blockC(iblock).tStimOneGratingContrast);  
        end
        if iblock == 1
            outS.trialsSinceReset = trs;
        else
            outS.trialsSinceReset = [outS.trialsSinceReset trs];
        end
    end
    for iF = 1:nF
        fN = cell2mat(fNames(iF,:));
        outS(1).(fN) = blockC(1).(fN);
        if size(blockC(1).(fN),2) == outS.trialsSinceReset(1);
            for iblock = 2:size(blockC,2)
                outS.(fN) = [outS.(fN) blockC(iblock).(fN)];
            end
        else
            outS.(fN) = blockC(1).(fN);
        end
    end