function tB = chkUniqueId(unqId)
% utility fn

% histed 121226

match = regexpi(unqId, '^[A-Za-z0-9_]*$');
if iscell(match)
    matchIx = ~cellfun('isempty', match);
else
    matchIx = match;
end
if ~matchIx
    error('UniqueId may only use numbers, characters, underscores (is %s)', ...
        unqId);
end

