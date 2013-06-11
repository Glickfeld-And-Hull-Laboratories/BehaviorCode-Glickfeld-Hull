function [outS ds0] = getBehavDataForDay(varargin)
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
    'RecData', [], ...
    'Subject', [], ...
    'DateStr', [], ...
    'WhichTrials', [], ...
    'DoCorrectEarlies', [], ...
    'DataIndex', 'max', ...
    'Debug', false, ...
    'MergeBlock1And2', false, ...
    'DoWhichIntensityIfBothPresent', '', ...  % empty means error; 'visual' or 'laser'
    'B2OneXTransformFH', [], ...
    'B2TwoXTransformFH', [], ...
    };

uo = stropt2struct(stropt_defaults(userDefs, varargin));

%% get data
lDebug = uo.Debug;

if isstruct(uo.Filename) 
    % backward compat
    assert(isempty(uo.Subject));
    ds = uo.Filename;
    ds = convertDataMatlab(ds);
elseif ~isempty(uo.RecData)
    assert(isfield(uo.RecData, 'MetaTags'), ...
        'Input does not look like a ''ns'' struct');
    ds = convertDataSerial(uo.RecData);
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
    ds =  mwLoadData(fName, dIndex, lDebug);
    if iscell(ds)
        ds = subConcatDataBlocks(ds);
    end
    ds = convertDataMatlab(ds);
end
if nargout > 1
    ds0 = ds;
end

%% backward compat on vectors
if ~isfield(ds, 'tLaserPowerMw') % before ~ Jan 2013
    ds.tLaserPowerMw = ds.laserPowerMw;
    ds.tGratingContrast = ds.gratingContrast;
    ds.tTotalReqHoldTimeMs = ds.reqHoldTimeMs;
end

%% process vectors in input  
laserPowerMw = ds.tLaserPowerMw;   % note these get converted to cell above in convertDataMatlab, I think MH 130402
gratingContrast = ds.tGratingContrast;
reqHoldTimesMs = ds.tTotalReqHoldTimeMs;

% handle case where one stimulus or the other was never presented
if isempty(laserPowerMw), laserPowerMw = nan([1 length(ds.trialOutcomeCell)]); end
if isempty(gratingContrast), gratingContrast = nan([1 length(ds.trialOutcomeCell)]); end

nLaserTrials = sum(~(isnan(laserPowerMw)|laserPowerMw==0));
nGratingTrials = sum(~(isnan(gratingContrast)|gratingContrast==0));
if nLaserTrials > 0 && nGratingTrials == 0
    intensityIsChr2 = true;
elseif nLaserTrials == 0 && nGratingTrials >0
    intensityIsChr2 = false;
else
    % both empty
    if isempty(uo.DoWhichIntensityIfBothPresent)
        error('Both both laser and visual stim present in data file and DoWhichIntensity is not set');
    elseif strcmpi(uo.DoWhichIntensityIfBothPresent, 'visual')
        intensityIsChr2 = false;
    elseif strcmpi(uo.DoWhichIntensityIfBothPresent, 'laser')
        intensityIsChr2 = true;
    else
        error('Unknown value for uo.DoWhichIntensityIfBothPresent');
    end
end
if intensityIsChr2 
    intensityV = laserPowerMw;
else
    intensityV = gratingContrast;
end

if isfield(ds, 'tBlock2TrialNumber') 
    block2V = ds.tBlock2TrialNumber;

    if uo.MergeBlock1And2 == true
        block2V = zeros(size(ds.tBlock2TrialNumber));
        warning('Overriding file block2 values and merging into one block');
    end
    
    if all(~ds.doBlock2) && (datenum(ds.startDateVec(1:3)) <= datenum([2012 06 07]))
        % kludge to deal with xml bug fixed on 120607 - tB2TN is set on all trials
        block2V = zeros(size(ds.tLaserPowerMw));
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
    intensityV = intensityV(desIx);

    block2V = block2V(desIx);
    reactTimesMs = ds.reactTimesMs(desIx);
    
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
assert(nBlock2Indices > 0 && nBlock2Indices <= 2);

intensityV = chop(intensityV,2);  % if you change the top value during session each power only accurate to 2 sig figs
successIx = strcmp(trialOutcomeCell, 'success');
earlyIx = strcmp(trialOutcomeCell, 'failure');
missedIx = strcmp(trialOutcomeCell, 'ignore');



%% consts 
allIntensities = unique(intensityV);
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
            xs = 1:max(ds.holdTimesMs+ds.reactTimeMs);
            earlyHoldV = ds.holdTimesMs(earlyIx);
            Fy = ksdensity(earlyHoldV, xs, 'function','pdf', 'width',300);
            
            nEarly = sum(earlyIx);
            nTrs = length(ds.tTotalReqHoldTimeMs);
            
            % for every correct find n trials and early rate in the previous rt window
            nC = sum(successIx);
            cList = find(successIx);
            for iC = 1:nC
                tTrN = cList(iC);
                tRt = ds.reactTimesMs(tTrN);
                tHold = ds.holdTimesMs(tTrN);
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
            
            estFalseCorrNPer100Ms = sum(ds.reactTimesMs >= -700 ...
                & ds.reactTimesMs <= tooFastMs) ...
                ./(700+tooFastMs)*100;
            estFalseCorrRatePer100Ms = estFalseCorrNPer100Ms ./ sum(nRawTot);
            
            ecCorrIx = ds.reactTimesMs >= tooFastMs & successIx;
            %reactWinLen = prctile(ds.reactTimesMs(corrIx), corrPrctile) - tooFastMs;
            
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

%% make cells
for iB=1:nBlock2Indices
    outS.percentsCorrectC{iB} = pctCorr(iB,1:nIntensities(iB));
end

%% output
fieldNames = { 'intensitiesC', 'percentsCorrect', 'intensityIsChr2', ...
    'nIntensities', 'nCorrPlusMissC', ...
    'successIx', 'earlyIx', 'missedIx', 'intensityV', ...
    'block2V', 'block2Indices', 'nBlock2Indices', ...
    'nMiss', 'nCorr', 'nCorrPlusMiss', 'whichTrials' };

nFs = length(fieldNames);
for iF=1:nFs
    tFN = fieldNames{iF};
    outS.(tFN) = eval([tFN ';']);
end
% backward compat - 
outS.reactTimesMs = ds.reactTimesMs;
% backward compat - before block2 implementation these were always non-cell
if nBlock2Indices == 1
    outS.intensV = intensitiesC{1};
    outS.pctCorr = percentsCorrect;
end

%% do any transforms 
fhList = { uo.B2OneXTransformFH, uo.B2TwoXTransformFH };
for iF=1:2
    tFH = fhList{iF};
    tB2N = iF;
        if ~isempty(tFH)
            assert(isa(tFH, 'function_handle'), ...
                'TransformFH must be a function handle');
            xPts = feval(tFH, outS.intensitiesC{tB2N}, tB2N-1, outS, ds);
            outS.intensitiesC{tB2N} = xPts;
            outS.intensitiyV = NaN*outS.intensityV; % if intensityV is needed edit this
        end
end







return
