function [ds, nDataBlocks nSelectedTrials] = mwLoadData(fName, dataIndex, debug)
%mwLoadData: Get data from matlab file saved by mworks
%
%  ds = mwLoadData(fName, dataN, debug)
%
%  ds is a struct if one block is returned, a cell vector containing structs if
%  more than one.
%
%  dataN can be 
%      numeric: return blocks with these indices
%      'last': return the last block
%      'max': return block with most trials
%      'max_N' return the N largest data blocks, keeping them in order
%      'all': return all
%  
% histed 110717

%120607 - removed trialThresholdPerChunk, added new dataN formats and
%    documented them.

%% arg processing
if nargin < 2 || isempty(dataIndex) || all(isnan(dataIndex)), dataIndex = 'last'; end
if nargin < 3 || isempty(debug); debug = false; end


%% get from disk
if ~exist(fName)
    error('MW:FileNotFound', 'missing data file %s', fName);
end
ds = load(fName);

%% figure out which chunk to load
if debug
    disp(sprintf('Filename %s', fName));
end

if isfield(ds, 'backup')
    ads = {ds.backup{:}, ds.input };
else
    ads = {ds.input};
end



    
nTrs = cellfun(@(x) length(x.holdStartsMs), ads);
nDataBlocks = length(ads);

if debug
    disp(sprintf('%d saved data blocks: nTrials %s', ...
        nDataBlocks, mat2str(nTrs)));
end

ds = {};
if isnumeric(dataIndex)
    nD = length(dataIndex);
    for iD = 1:nD
        ds{iD} = ads{dataIndex(iD)};
    end
elseif ischar(dataIndex)
    [maxMatchSt,~,~,~,toks] = regexpi(dataIndex, 'max_?([0-9]*)');
    if strcmpi(dataIndex, 'last')
        ds = ads(end);  
    elseif strcmpi(dataIndex, 'all')
        ds = ads;
    elseif ~isempty(maxMatchSt)
        [sortVals sortNs] = sort(nTrs);
        sortVals = fliplr(sortVals);
        sortNs = fliplr(sortNs);

        % choose max indices
        tN = str2double(toks{1}{1});
        if isnan(tN), tN = 1; end  % 'max' case
        maxNs = 1:tN;
        desNs = sort(sortNs(maxNs));
        ds = ads(desNs);
       
        if debug
            disp(sprintf('Using max trials: chunks %s', mat2str(desNs)));
        end
    end
end
if isempty(ds)
    error('invalid dataIndex: %s', mat2str(dataIndex));
end

desTrNs = cellfun(@(x) length(x.trialOutcomeCell), ds);
if debug
    disp(sprintf('DataIndex %s: Selected chunk(s): %s trials', ...
        mat2str(dataIndex), mat2str(desTrNs)));
end

%% format output
if nargout >2
    nSelectedTrials = desTrNs;
end
if iscell(ds) && length(ds) == 1
    ds = ds{1}; % unpack cell
end 

