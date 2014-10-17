function [outS ds0] = getBehavDataForDayLG(varargin)
% encapsulate code to get psychometric data for a given day
%
% If only one block2 level is found in file, this 
%  userDefs = { ...
%    'Filename', [], ...
%    'Subject', [], ...
%    'DateStr', [], ...
%    'WhichTrials', [], ...
%    'DoCorrectEarlies', [], ...
%    'DataIndex', 'max' };


userDefs = { ...
    'Filename', [], ...
    'Subject', [], ...
    'DateStr', [], ...
    'WhichTrials', [], ...
    'DoCorrectEarlies', [], ...
    'DataIndex', 'max', ...
    'Debug', false, ...
    'MergeBlock1And2', false, ...
    'DoBlock1Only', false, ...
    'SplitBlock1', false, ...
    'MergeMats', false,};

uo = stropt2struct(stropt_defaults(userDefs, varargin));

%% get data
lDebug = uo.Debug;

if isstruct(uo.Filename)
    assert(isempty(uo.Subject));
    ds = uo.Filename;
else
    % read from disk

    % construct fname from subjnum
    if isempty(uo.Filename)
        rtc = behavConstsHADC8;
        fName = rtc.computeFName(uo.Subject, uo.DateStr);
    else
        fName = uo.Filename;
    end

    if ischar(uo.DataIndex)
        if strcmpi(uo.DataIndex, 'max')
            dIndex = 'max';
        else
            try
                dIndex = eval(uo.DataIndex);  
            catch mexc
                error('Invalid DataIndex field, in converting to mat: %s', uo.DataIndex);
                rethrow mexc;
            end
        end
    else
        dIndex = uo.DataIndex;
    end
    %compatibility with new naming HHMM
    if ~exist(fName)
        rtc = behavConstsHADC8;
        n  = dir([fName(1:45) '*']);
        if size(n,1) == 1
            fName = fullfile(rtc.pathStr, n.name);
            ds =  mwLoadData(fName, dIndex, lDebug);
        elseif size(n,1) > 1
            if uo.MergeMats == 1
                for ifile = 1:size(n,1)
                    fName = fullfile(rtc.pathStr, n(ifile).name);
                    ds(ifile) =  mwLoadData(fName, dIndex, lDebug);
                    ds = concatenateDataBlocks(ds);
                end
            else
                error('Too many mat files- need to choose');
            end              
        end
    else
         ds =  mwLoadData(fName, dIndex, lDebug);
    end
    if iscell(ds)
        ds = subConcatenateDataBlocks(ds);
    end
end

%backwards compatibility
if isfield(ds, 'stimOnTimeMs')
    ds.gratingDurationMs = ds.stimOnTimeMs;
end

if nargout > 1
    ds0 = ds;
end

        

%% process vectors in input
nTrTemp = size(ds.trialOutcomeCell);
vFields = { 'tLaserPowerMw', 'tGratingContrast', 'tBaseGratingContrast', 'tTotalReqHoldTimeMs', ...
    'tGratingDirectionDeg', ...
    'tTotalReqHoldTimeMs', 'holdTimesMs' };
for iV=1:length(vFields)
    if isfield(ds, vFields{iV})
        eval(sprintf('%s = celleqel2mat_padded(ds.%s);', vFields{iV}, vFields{iV}));
    else
        eval(sprintf('%s = repmat(NaN, [1 nTrTemp]);', vFields{iV}));
    end
end

laserPowerMw = tLaserPowerMw;
gratingContrast = tGratingContrast;
gratingDirectionDeg = tGratingDirectionDeg;
reqHoldTimesMs = tTotalReqHoldTimeMs;
baseGratingContrast = tBaseGratingContrast;

nLaserTrials = sum(~(isnan(laserPowerMw)|laserPowerMw==0));
nContrastTrials = sum(find(gratingContrast>min(tGratingContrast,[],2)));
nOriTrials = sum(find(gratingDirectionDeg>min(tGratingDirectionDeg,[],2)));
nGratingTrials = nContrastTrials + nOriTrials;

if nLaserTrials > 0 && nGratingTrials == 0
    intensityIsChr2 = true;
end
if nLaserTrials == 0 && nGratingTrials >0
    intensityIsChr2 = false;
    if nOriTrials>0
        doContrast = false;
        doOri = true;
    else
        doContrast = true;
        doOri = false;
    end
elseif ~isfield(ds, 'tGratingContrast')
    %hack for changed names
    gratingContrast = celleqel2mat_padded(ds.gratingContrast);
    nContrastTrials = sum(find(gratingContrast>min(gratingContrast,[],2)));
    nOriTrials = sum(find(tGratingDirectionDeg>min(tGratingDirectionDeg,[],2)));
    nGratingTrials = nContrastTrials + nOriTrials;
    if nGratingTrials > 0
        intensityIsChr2 = false;
        if nOriTrials>0
            doContrast = false;
            doOri = true;
        else
            doContrast = true;
            doOri = false;
        end
    end
    elseif isfield(ds, 'tGratingContrast')
    % both empty
    error('Cannot handle such data files w/ both laser and grating yet - check code');
end

if intensityIsChr2 
    intensityV = laserPowerMw;
elseif nOriTrials>0
    intensityV = gratingDirectionDeg;
elseif nContrastTrials>0
    intensityV = gratingContrast;
end

reactTimesMs = double(cellvect2mat_padded(ds.reactTimesMs));
if isfield(ds, 'tBlock2TrialNumber') 
    block2V = cellvect2mat_padded(ds.tBlock2TrialNumber);

    if uo.MergeBlock1And2 == true
        block2V = zeros(size(ds.tBlock2TrialNumber));
        warning('Overriding file block2 values and merging into one block');
    end
    
    if ~ds.doBlock2 && (datenum(ds.startDateVec(1:3)) <= datenum([2012 06 07]))
        % kludge to deal with xml bug fixed on 120607 - tB2TN is set on all trials
        block2V = zeros(size(laserPowerMw));
    end
    
    if uo.SplitBlock1 == true
        block2V = zeros(size(ds.tBlock2TrialNumber));
        trialLaserOn = ds.trialLaserOnTimeMs;
        trialLaserOff = ds.trialLaserOffTimeMs;
        block2V(find(rem(reqHoldTimesMs,(trialLaserOn+trialLaserOff))<=150)) = 1;
        block2V(find(rem(reqHoldTimesMs,(trialLaserOn+trialLaserOff))>150)) = 0;
        block2V(find(rem(reqHoldTimesMs,(trialLaserOn+trialLaserOff))>350)) = 1;
    end
else
    block2V = zeros(size(intensityV));  % no values in file, create a single block
end
trialOutcomeCell = ds.trialOutcomeCell;  % copy so we can trim if necessary below
%% trim trials if requested
if ~isempty(uo.WhichTrials)
    if ischar(uo.WhichTrials)
        try
            whichTrials = eval(uo.WhichTrials);
        catch mexc
            warning('Error parsing text WhichTrials specification: %s', uo.WhichTrials);
            rethrow(mexc);
        end
    else
        whichTrials = uo.WhichTrials;
    end
    
    desIx = false(size(trialOutcomeCell));
    desIx(whichTrials) = true;

    try  % use direct indexing to figure out if any indices are invalid
        trialOutcomeCell = trialOutcomeCell(desIx);  
    catch
        nTrials = length(trialOutcomeCell);
        if max(uo.WhichTrials) > nTrials
            error('Asked for a trial beyond nTrials: max asked %d, nTrials %d', ...
                max(uo.WhichTrials), nTrials);
        else
            error('Bad WhichTrials specification: %s', mat2str(uo.WhichTrials));
        end
    end
    if length(trialOutcomeCell) == 0
        error('No trials remaining after selection');
    end

    
    laserPowerMw = laserPowerMw(desIx);
    gratingContrast = gratingContrast(desIx);
    gratingDirectionDeg = gratingDirectionDeg(desIx);
    intensityV = intensityV(desIx);

    block2V = block2V(desIx);
    reactTimesMs = reactTimesMs(desIx);
    
    nDes = sum(desIx);
    nTotal = length(desIx);
    fprintf(1,'%s: trimmed %d trials; requested via WhichTrials (%3.1f%%, %d to %d trs)\n', ...
            mfilename, ...
            nTotal-nDes, (nTotal-nDes)/nTotal*100, nTotal, nDes);
    
else
    desIx = true(size(trialOutcomeCell));
    whichTrials = find(desIx);
end

%% process vectors post-trim
if all(isnan(block2V))
    block2V = zeros(size(block2V));
end
block2Indices = unique(block2V(~isnan(block2V)));
nBlock2Indices = length(block2Indices);
    
if uo.DoBlock1Only == true
    nBlock2Indices = 1;  % no values in file, create a single block
end

assert(nBlock2Indices > 0 && nBlock2Indices <= 2);

intensityV = chop(intensityV,2);  % if you change the top value during session each power only accurate to 2 sig figs
successIx = strcmp(trialOutcomeCell, 'success');
earlyIx = strcmp(trialOutcomeCell, 'failure');
missedIx = strcmp(trialOutcomeCell, 'ignore');



%% consts 
for iB = 1:nBlock2Indices
    tBN = block2Indices(iB);
    b2Ix = block2V == tBN;
    
    intensitiesC{iB} = unique(intensityV(b2Ix));
    nIntensities(iB) = length(intensitiesC{iB});
    
    intensityNums = block2V*NaN;
    %% compute number of corrects and rts
    for iI = 1:nIntensities(iB)
        tI = intensitiesC{iB}(iI);
        iIx = intensityV == tI;
        
        nCorr(iB, iI) = sum(iIx & successIx & b2Ix);
        nEarly(iB, iI) = sum(iIx & earlyIx & b2Ix);
        nMiss(iB, iI) = sum(iIx & missedIx & b2Ix);
        nRawTot(iB, iI) = sum(iIx & b2Ix);

        tV = reactTimesMs(iIx & successIx & b2Ix);
        if length(tV) == 0, 
            assert( sum(iIx & b2Ix) > 0 );
            assert( sum(iIx & b2Ix & successIx) == 0);
            tV = NaN; 
        end 
        outS.reactTimesByPower{iI,iB} = tV;
        outS.reactTimeMean(iI,iB) = mean(tV);
        outS.reactTimeStd(iI,iB) = std(tV);
        outS.reactTimeSEM(iI,iB) = std(tV) ./ sqrt(length(tV));
        
        outS.intensityNums(iIx) = iI;
        
    end
    nCorrPlusMissC{iB} = nCorr(iB,1:nIntensities(iB)) + nMiss(iB,1:nIntensities(iB));
end
nCorrPlusMiss = nCorr+nMiss;





%% correct for early rate 
if uo.DoCorrectEarlies

    whichEarlyMethod = 'earlyPDF'; %'earlyPDF';  % or 'old'
    
    switch whichEarlyMethod
        case 'earlyPDF'
            %% method 1 - based on estimating early rate and extrap over actual RT window
            realTooFast = 125; % min RT to any stimulus
            % first find early PDF
            xs = 1:max(holdTimesMs+ds.reactTimeMs);
            earlyHoldV = holdTimesMs(earlyIx);
            Fy = ksdensity(earlyHoldV, xs, 'function','pdf', 'width',300);
            
            nEarly = sum(earlyIx);
            nTrs = length(reqHoldTimesMs);
            
            % for every correct find n trials and early rate in the previous rt window
            nC = sum(successIx);
            cList = find(successIx);
            for iC = 1:nC
                tTrN = cList(iC);
                tRt = reactTimesMs(tTrN);
                tHold = holdTimesMs(tTrN);
                earlyCdfIx = xs >= tHold-tRt-realTooFast & xs < tHold;  % probably don't need the toofast addition but it makes us more likely to reject trials/conservative
                earlyRateInRtWin = sum(Fy(earlyCdfIx));
                probEarly(iC) = earlyRateInRtWin;
            end
            probEarlyOrNone = zeros(size(successIx));
            probEarlyOrNone(successIx) = probEarly;
            
            for iB=1:nBlock2Indices
                for iI=1:nIntensities(iB)
                    tI = intensitiesC{iB}(iI);
                    iIx = intensityV == tI;
                    b2Ix = block2V == block2Indices(iB);
                    estFalseCorrPerImage(iB,iI) = sum(probEarlyOrNone(iIx & b2Ix));%
                end
            end
                
        case 'old'
            %% old method
            tooFastMs = 125; % use a 'true' toofasttime here % ds.tooFastTimeMs;
            corrPrctile = 90;
            
            estFalseCorrNPer100Ms = sum(reactTimesMs >= -700 ...
                & reactTimesMs <= tooFastMs) ...
                ./(700+tooFastMs)*100;
            estFalseCorrRatePer100Ms = estFalseCorrNPer100Ms ./ sum(nRawTot);
            
            ecCorrIx = reactTimesMs >= tooFastMs & successIx;
            %reactWinLen = prctile(reactTimesMs(corrIx), corrPrctile) - tooFastMs;
            
            reactWinLen = ds.reactTimeMs - tooFastMs;
            if length(reactWinLen) > 1
                assert(range(reactWinLen) < 0.1*max(reactWinLen), 'bug- RT window changed a lot?');
                reactWinLen = min(reactWinLen); % err on side of underestimation
            end
            
            estFalseCorrPerImage = nCorrPlusMiss.*estFalseCorrRatePer100Ms.*reactWinLen/100;
    end
            
    nCorrOrig = nCorr;
    nCPMOrig = nCorrPlusMiss;
    nCorr = round(nCorr - estFalseCorrPerImage);
    nCorrPlusMiss = round(nCorrPlusMiss - estFalseCorrPerImage);
    
    % set output
    outS.estFalseCorrPerImage = fliplr(estFalseCorrPerImage);
    outS.nCorrOrig = nCorrOrig;
    
    if any(nCorr < 0)
        nCorr(nCorr<0) = 0;
    end
  
    disp(sprintf('Removed total %d estimated false corrects (method %s)- by image %s', ...
        sum(round(outS.estFalseCorrPerImage(:))), whichEarlyMethod, mat2str(chop(outS.estFalseCorrPerImage,2))));
else
    outS.estFalseCorrPerImage = [];
end


pctCorr = nCorr./nCorrPlusMiss;

pctCorr(pctCorr==1) = pctCorr(pctCorr==1)-10*eps;

percentsCorrect = pctCorr;
%% correct intensities
if doContrast
    if max(baseGratingContrast,[],2)>0
        if ds.gratingMaxContrastStep<0
            for iB = 1:nBlock2Indices
                intensitiesC{iB} = max(baseGratingContrast,[],2)-intensitiesC{iB};
            end
        else
            for iB = 1:nBlock2Indices
                intensitiesC{iB} = intensitiesC{iB}-max(baseGratingContrast,[],2);
            end
        end
    end
end
    
%% make cells
outS.percentsCorrectC{1} = pctCorr(1,~isnan(pctCorr(1,:)));
if nBlock2Indices == 2
    outS.percentsCorrectC{2} = pctCorr(2,~isnan(pctCorr(2,:)));
end

%% output
fieldNames = { 'intensitiesC', 'percentsCorrect', 'intensityIsChr2', ...
    'nIntensities', 'nCorrPlusMissC', ...
    'doContrast', 'doOri', ...
    'successIx', 'earlyIx', 'missedIx', 'intensityV', 'reactTimesMs', ...
    'block2V', 'block2Indices', 'nBlock2Indices', ...
    'nMiss', 'nCorr', 'nCorrPlusMiss', 'whichTrials' };

nFs = length(fieldNames);
for iF=1:nFs
    tFN = fieldNames{iF};
    outS.(tFN) = eval([tFN ';']);
end
% backward compat - before block2 implementation these were always non-cell
if nBlock2Indices == 1
    outS.intensV = intensitiesC{1};
    outS.pctCorr = percentsCorrect;
end



function outS = subConcatenateDataBlocks(blockC)

outS = struct([]);

% going to be pretty hackish
fNames = fieldnames(blockC{1});
nF = length(fNames);
for iF = 1:nF
    tDs1 = blockC{1};
    tFName = fNames{iF};
    tFV1 = tDs1.(tFName);
    
    lenList = cellfun(@(x) length(x.(tFName)), blockC);
    
    newC = cellfun(@(x) x.(tFName), blockC, 'UniformOutput', false);

    if isnumeric(tFV1) || islogical(tFV1)
        newV = cat(2, newC{:});
                
        % reduce to one element if all the same
        if all(lenList) == 1
            if all(newV == newV(1))
                newV = newV(1);
            end
        end
        
        outS(1).(tFName) = newV;
    elseif ischar(tFV1)
        if all(strcmp(newC, newC{1}));
            newC = newC{1};
        end
        outS(1).(tFName) = newC;
    elseif iscell(tFV1)
        % check same
        if all(cellfun(@(x) isequalwithequalnans(x.(tFName), tFV1), blockC))
            newC = newC{1};
        elseif isvector(tFV1)
            newC = cat(2, newC{:});
        % else use newC as is
        end
        outS(1).(tFName) = newC;
    elseif isstruct(tFV1)
        assert(strcmp(tFName, 'firstTrConsts'));
        outS(1).(tFName) = blockC{1}.(tFName);
    else
        error('missing field');
    end
end


%% combining multiple mat files

function outS = concatenateDataBlocks(blockC)
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

