function bxStruct = bxFrmXLS(mouseName,bxDataInfo)
    rc = behavConstsAV;
    nexp = size(bxDataInfo,1);
    bxStruct = struct;
    bxParams_FSAV;
    for iexp = 1:nexp
        fprintf(num2str(iexp))
        expDate = bxDataInfo(iexp).DateStr;
        if isa(expDate,'double')
            expDate = num2str(expDate);
        end
        runs = bxDataInfo(iexp).ChooseMatFile;
        if isempty(bxDataInfo(iexp).TrialRangeToUse)
            trialRange = [];
        else
            trialRange = eval(bxDataInfo(iexp).TrialRangeToUse);
        end

        mworksFilesDir = dir(fullfile(rc.behavData,...
            sprintf('data-i%s-%s-*',mouseName,expDate)));
        if isempty(runs) | isnan(runs)
            runs = nan;
        elseif isa(runs,'double')
            runs = runs;
        else
            runs = eval(runs);
        end
        mworks = [];
        if isnan(runs)
            mworksTime = mworksFilesDir(1).name(end-7:end-4);
            mworks = loadMworksFile(mouseName,expDate,mworksTime);
        else
            mworksFilesDir = mworksFilesDir(runs);
            for irun = 1:length(runs)
                mworksTime = mworksFilesDir(1).name(end-7:end-4);
                mworks_temp = loadMworksFile(mouseName,expDate,mworksTime);
                mworks = [mworks,mworks_temp];
                clear mworks_temp
            end
            mworks = concatenateDataBlocks(mworks);
        end

        if ~isempty(trialRange)
            if trialRange(end) > length(mworks.trialOutcomeCell)
                lastTrial = find(trialRange > length(mworks.trialOutcomeCell),1)-1;
                trialRange = trialRange(1:lastTrial);
            end
            trialRangeInd = false(1,length(mworks.trialOutcomeCell));
            trialRangeInd(trialRange) = true;
            mworks = trialChopper(mworks,trialRangeInd);
        end
        
        
        if isfield(mworks,'nFramesOn')
            frameRate = mworks.frameRateHz;
            
            trialStartFr = celleqel2mat_padded(mworks.cFirstStim);
            trialEndFr = celleqel2mat_padded(mworks.cLeverUp);
            trialTimeFr = trialEndFr - trialStartFr;
            
            trialTimeMs = (trialTimeFr./frameRate).*1000;
        else
            leverUpTime = cell2mat(mworks.leverUpTimeMs);
            leverDownTime = cell2mat(mworks.leverDownTimeMs);
            trialTimeMs = double(leverUpTime - leverDownTime);
        end
        
        tooShortTrials = trialTimeMs < minTrialLengthMs;
        if any(tooShortTrials)
            mworks = trialChopper(mworks,~tooShortTrials);
            if isfield(mworks,'nFramesOn')
                frameRate = mworks.frameRateHz;

                trialStartFr = celleqel2mat_padded(mworks.cFirstStim);
                trialEndFr = celleqel2mat_padded(mworks.cLeverUp);
                trialTimeFr = trialEndFr - trialStartFr;

                trialTimeMs = (trialTimeFr./frameRate).*1000;
            else
                leverUpTime = cell2mat(mworks.leverUpTimeMs);
                leverDownTime = cell2mat(mworks.leverDownTimeMs);
                trialTimeMs = double(leverUpTime - leverDownTime);
            end
        end
        
        nCyc = double(cell2mat(mworks.nCyclesOn));
        tOn = double(mworks.stimOnTimeMs);
        tOff = double(mworks.stimOffTimeMs);

        valTargetOnTime =(tOn+tOff)*nCyc;
        valReactTimeCalc = trialTimeMs - valTargetOnTime;
        
        nTrials = length(mworks.trialOutcomeCell);
        fa = strcmp(mworks.trialOutcomeCell,'failure');
        miss = strcmp(mworks.trialOutcomeCell,'ignore');
        hit = strcmp(mworks.trialOutcomeCell,'success');
        
%         fidgets = trialTimeMs < (tOn+tOff+100);
%         fa(fidgets) = false; 
        
        pctEarly = sum(fa)./nTrials;
        
        tVisTargets = round(double(cell2mat_padded(mworks.tGratingDirectionDeg)),2,'significant');
        visTargets = unique(tVisTargets);
        visTargets = visTargets(2:end);        
        visHR = zeros(1,length(visTargets));
        nVis = zeros(1,length(visTargets));
        for i = 1:length(visTargets)
            ind = tVisTargets == visTargets(i);
            visHR(i) = sum(ind & hit) ./ sum(ind & (hit | miss));
            nVis(i) = sum(ind & (hit | miss));
        end
        
        tAudTargets = round(double(celleqel2mat_padded(...
            mworks.tSoundTargetAmplitude)),2,'significant');
        audTargets = unique(tAudTargets);
        audTargets = audTargets(2:end);
        audHR = zeros(1,length(audTargets));
        nAud = zeros(1,length(audTargets));
        for i = 1:length(audTargets)
            ind = tAudTargets == audTargets(i);
            audHR(i) = sum(ind & hit) ./ sum(ind & (hit | miss));
            nAud(i) = sum(ind & (hit | miss));
        end
        
        bxStruct(iexp).date = expDate;
        bxStruct(iexp).sn = mouseName;
        bxStruct(iexp).pctEarly = pctEarly;
        bxStruct(iexp).visHR = visHR;
        bxStruct(iexp).audHR = audHR;
        bxStruct(iexp).hit = hit;
        bxStruct(iexp).miss = miss;
        bxStruct(iexp).fa = fa;
        bxStruct(iexp).tVisTargets = tVisTargets;
        bxStruct(iexp).tAudTargets = tAudTargets;
        
        bxStruct(iexp).stimOnTime = mworks.stimOnTimeMs;
        bxStruct(iexp).stimOffTime = mworks.stimOffTimeMs;
        bxStruct(iexp).trLengthMs = trialTimeMs;
        bxStruct(iexp).nCycles = nCyc;
        bxStruct(iexp).valReact = valReactTimeCalc;
        bxStruct(iexp).size = mworks.gratingHeightDeg;
        
        if mworks.doShortCatchTrial
            bxStruct(iexp).tInvVisTargets = cell2mat_padded(...
                mworks.tCatchGratingDirectionDeg);
            bxStruct(iexp).tInvAudTargets = celleqel2mat_padded(...
                mworks.tSoundCatchAmplitude);
            if isempty(bxStruct(iexp).tInvVisTargets)
                bxStruct(iexp).tInvVisTargets = nan(1,length(bxStruct(iexp).tInvAudTargets));
            elseif isempty(bxStruct(iexp).tInvAudTargets)
                bxStruct(iexp).tInvAudTargets = nan(1,length(bxStruct(iexp).tInvVisTargets));
            end
            if mworks.doShortCatchReward == 0
                disp(['Catch not rewared - ' expDate '-' mouseName])
            end
            
            if isfield(mworks,'cCatchOn')
                cCyc = celleqel2mat_padded(mworks.cCatchOn);
                
                catchTrialFr = cCyc - trialStartFr;
                catchReactFr = trialTimeFr - catchTrialFr;
                
                tooFastFr = round(frameRate*((mworks.tooFastTimeMs)/1000));
                reactTimeFr = round(frameRate*((mworks.reactTimeMs)/1000));
                
                catchOutcome = cell(1,nTrials);
                outInd = catchReactFr > tooFastFr & catchReactFr < reactTimeFr;
                catchOutcome(outInd) = {'FA'};
                outInd = catchReactFr > reactTimeFr;
                catchOutcome(outInd) = {'CR'};
                outInd = catchReactFr < (0 +tooFastFr);
                catchOutcome(outInd) = {'failure'};
                bxStruct(iexp).invHit = strcmp(catchOutcome,'FA');
                bxStruct(iexp).invMiss = strcmp(catchOutcome,'CR');
                
                if iscell(mworks.nFramesOn)
                    frOn = unique(cell2mat(mworks.nFramesOn));
                    frOff = unique(cell2mat(mworks.nFramesOff));
                else
                    frOn = mworks.nFramesOn;
                    frOff = mworks.nFramesOff;
                end
                catchTimeMs = (catchTrialFr./frameRate)*1000;
                catchCycCalc = double(round(catchTrialFr./(frOn+frOff)));
                invReactTimeCalc = trialTimeMs - catchTimeMs;
                
            else
                bxStruct(iexp).invHit = strcmp(mworks.catchTrialOutcomeCell,'FA');
                bxStruct(iexp).invMiss = strcmp(mworks.catchTrialOutcomeCell,'CR');
                catchTimeMs = celleqel2mat_padded(mworks.tCatchTimeMs) - ...
                    double(leverDownTime);
                catchCycCalc = double(round(catchTimeMs./(tOn+tOff)));
%                 invTargetOnTime = (tOn+tOff)*catchCycCalc;
%                 invReactTimeCalc = trialTimeMs - invTargetOnTime;
                invReactTimeCalc = trialTimeMs - catchTimeMs;
            end
            if any(bxStruct(iexp).tInvVisTargets > 0)
                if any(bxStruct(iexp).tInvAudTargets > 0)
                    bxStruct(iexp).invType = 'both';
                else
                    bxStruct(iexp).invType = 'vis';
                end
            elseif any(bxStruct(iexp).tInvAudTargets > 0)
                bxStruct(iexp).invType = 'aud';
            end
            
            bxStruct(iexp).catchCyc = catchCycCalc;
            bxStruct(iexp).invReact = invReactTimeCalc;
        else
            bxStruct(iexp).invType = nan;
            bxStruct(iexp).invHit = nan(1,nTrials);
            bxStruct(iexp).invMiss = nan(1,nTrials);
            bxStruct(iexp).catchCyc = nan(1,nTrials);
            bxStruct(iexp).invReact = nan(1,nTrials);
            
        end 
            
    end
end
