function bxStruct = createBxStruct(exptStruct)
    rc = behavConstsAV;
    nexp = size(exptStruct,2);
    bxStruct = struct;
    bxParams_FSAV_attnV1ms;
    for iexp = 1:nexp
        fprintf(num2str(iexp))
        mouseName = exptStruct(iexp).SubNum;
        expDate = exptStruct(iexp).date;
        runs = exptStruct(iexp).runs;
        if isnan(exptStruct(iexp).trial_range)
            trialRange = [];
        else
            trialRange = exptStruct(iexp).trial_range;
        end
        nrun = size(runs,1);
        mworks = [];
        for irun = 1:nrun
            mw_temp = loadMworksFile(mouseName,expDate,exptStruct(iexp).time_mat(irun,:));
            if irun > 1 & ~isempty(setdiff(fieldnames(mworks),fieldnames(mw_temp)))
                extra_fn = setdiff(fieldnames(mworks),fieldnames(mw_temp));
                for ifn = 1:length(extra_fn)
                    if length(eval(['mworks.' extra_fn{ifn}])) > 1
                        mw_temp.(extra_fn{ifn}) = cell(1,mw_temp.trialSinceReset);
                    else
                        mw_temp.(extra_fn{ifn}) = nan;
                    end
                end
            else
                mworks = [mworks,mw_temp];
            end
        end
        mworks = concatenateDataBlocks(mworks);

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

        earlyHits_vis = hit & valReactTimeCalc < visRTwindow(1) & tVisTargets > 0;
        lateHits_vis = hit & valReactTimeCalc > visRTwindow(2) & tVisTargets > 0;
        
        hit(earlyHits_vis | lateHits_vis) = false;
        miss(lateHits_vis) = true;
        fa(earlyHits_vis) = true;
        
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
        

        earlyHits_aud = hit & valReactTimeCalc < visRTwindow(1) & tAudTargets > 0;
        lateHits_aud = hit & valReactTimeCalc > visRTwindow(2) & tAudTargets > 0;
        
        hit(earlyHits_aud | lateHits_vis) = false;
        miss(lateHits_aud) = true;
        fa(earlyHits_aud) = true;

        if any(earlyHits_vis | lateHits_vis |earlyHits_aud | lateHits_aud)
            fprintf('Expt %s: %s/%s trials changed due to RT cutoff\n',num2str(iexp),...
                num2str(sum(earlyHits_vis | lateHits_vis |earlyHits_aud | lateHits_aud)),...
                num2str(length(hit)))
        end

        bxStruct(iexp).date = expDate;
        bxStruct(iexp).sn = mouseName;
        bxStruct(iexp).pctEarly = pctEarly;
        bxStruct(iexp).trialsPerRun = mworks.trialsSinceReset;
        bxStruct(iexp).isTrueVisTrial = cell2mat(mworks.tBlock2TrialNumber) == 0;
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

            earlyHit_invVis = bxStruct(iexp).invHit ...
                & bxStruct(iexp).tInvVisTargets > 0 ...
                & invReactTimeCalc < visRTwindow(1);
            lateHit_invVis = bxStruct(iexp).invHit ...
                & bxStruct(iexp).tInvVisTargets > 0 ...
                & invReactTimeCalc > visRTwindow(2);
            earlyHit_invAud = bxStruct(iexp).invHit ...
                & bxStruct(iexp).tInvAudTargets > 0 ...
                & invReactTimeCalc < visRTwindow(1);
            lateHit_invAud = bxStruct(iexp).invHit ...
                & bxStruct(iexp).tInvAudTargets > 0 ...
                & invReactTimeCalc > visRTwindow(2);

            bxStruct(iexp).invHit(earlyHit_invVis | lateHit_invVis ...
                | earlyHit_invAud | lateHit_invAud) = false;
            bxStruct(iexp).invMiss(lateHit_invVis | lateHit_invAud) = true;

            
            if any(earlyHit_invVis | lateHit_invVis |earlyHit_invAud | lateHit_invAud)
                fprintf('Expt %s: %s/%s invalid trials changed due to RT cutoff\n',num2str(iexp),...
                    num2str(sum(earlyHit_invVis | lateHit_invVis |earlyHit_invAud | lateHit_invAud)),...
                    num2str(sum(bxStruct(iexp).tInvVisTargets > 0 |  bxStruct(iexp).tInvAudTargets > 0)))
            end
        else
            bxStruct(iexp).invType = nan;
            bxStruct(iexp).invHit = nan(1,nTrials);
            bxStruct(iexp).invMiss = nan(1,nTrials);
            bxStruct(iexp).catchCyc = nan(1,nTrials);
            bxStruct(iexp).invReact = nan(1,nTrials);
            
        end 
            
    end
end
