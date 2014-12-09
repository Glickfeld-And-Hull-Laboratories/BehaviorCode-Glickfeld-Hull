function [outStruct] = trialChopper(input, tr, randTrials)
%% TrialChopper: [outStruct] = trialChopper(input, [tr], randTrials)
% input =  data structure to be indexed
% tr = [trialStart trialEnd], or a boolean index
% randTrials = number of trials to be selected between trialStart and
% trialEnd (optional)
%
% Designed to index a range of trials within a data structure, then randomly select
% randTrials number of trials from that range, then output the reformatted
% structure.
if nargin>1,
    if length(tr)==2 && issorted(tr),
        trStart = tr(1);
        trEnd = tr(2);
        iX = trStart:trEnd;
    else
        iX = logical(tr);
    end
    if nargin>2,
        rTs = randTrials;
        rg = trEnd-trStart;
        randIx = randperm(rg, rTs);
    else
        rTs = [];
    end
end
            
ds = input;
outputStruct = ds;
fieldNames = fieldnames(ds);
nFNs = length(fieldNames);
nTrs = length(ds.trialOutcomeCell);

for i=1:nFNs,
    a = strcat('input.', fieldNames(i));
    field = eval(a{1});
    if length(field)==nTrs && ~strcmp(fieldNames(i), 'savedEvents'),
        field = field(iX);
        outStruct.(fieldNames{i}) = field;
        if ~isempty(rTs),
            outStruct.(fieldNames{i}) = field(randIx);
        end
    else
        outStruct.(fieldNames{i}) = field;
    end    
end
outStruct.trialSinceReset = 1;