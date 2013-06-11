function endStateS = mwGetEndState(ds)
%mwGetEndState: unpack events from matlab file and return final state each var
%
% note: 121203- not sure this works well but it may be a useful example so
% I'm leaving it here.
%
% histed 110717

if ~isfield(ds, 'eventCodecs') || isempty(ds.eventCodecs) || isempty(ds.savedEvents{1})
    fprintf(1, '%s: No event data found in file, returning empty\n', mfilename);
    endStateS = [];
    return
end

lCodec = ds.eventCodecs{end};


%% testing
b = cellfun(@struct2cell, ds.savedEvents, 'UniformOutput', false);
b1 = cat(2,b{:}); % rows: event_code, time_us, data

codeNums = cat(2,b1{1,:});
%codeTimes = cat(2,b1{2,:});
%codeDataC = b1{3,:});  % draw from orig cell

allCodeNames = {lCodec.tagname};
nCodes = length(allCodeNames);
for iC = 1:nCodes
    tName = lCodec(iC).tagname;
    if tName(1) == '#';
        tName = cat(2,'x_', tName(2:end)); % prepend x for struct field name
    end
    
    tVal = lCodec(iC).code;
    
    desIx = codeNums == tVal;
    desN = find(desIx, 1, 'last');
    
    if isempty(desN)
        tData = [];
    else
        tData = b1{3,desN};
    end

    endStateS(1).(tName) = tData;
end

    
% $$$ %% old code, slow!!
% $$$ lCCodes = cat(2, lCodec.code);
% $$$ lCNames = {lCodec.tagname};
% $$$ [codeToStr{lCCodes}] = deal(lCNames{:});
% $$$ 
% $$$ nE = length(ds.savedEvents);
% $$$ endStateS = struct([]);
% $$$ 
% $$$ 
% $$$ keyboard
% $$$ for iE = 1:nE
% $$$     tE = ds.savedEvents{iE};
% $$$     tCodes = cat(2,tE.event_code);
% $$$     tNames = codeToStr(tCodes);
% $$$     nN = length(tNames);
% $$$     for iN = 1:nN
% $$$         tN = tNames{iN};
% $$$         endStateS(1).(tN) = tE(iN).data;
% $$$     end
% $$$ end
% $$$ 
