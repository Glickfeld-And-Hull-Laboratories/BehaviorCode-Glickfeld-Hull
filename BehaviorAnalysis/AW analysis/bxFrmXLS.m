function bxStruct = bxFrmXLS(mouseName,bxDataInfo)
    rc = behavConstsAV;
    nexp = size(bxDataInfo,1);
    bxStruct = struct;
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
                trialRange(end) = length(mworks.trialOutcomeCell);
            end
            mworks = trialChopper(mworks,[trialRange(1) trialRange(end)]);
        end
        
        nTrials = length(mworks.trialOutcomeCell);
        fa = strcmp(mworks.trialOutcomeCell,'failure');
        miss = strcmp(mworks.trialOutcomeCell,'ignore');
        hit = strcmp(mworks.trialOutcomeCell,'success');
        
        pctEarly = sum(fa)./nTrials;
        
        tVisTargets = cell2mat_padded(mworks.tGratingDirectionDeg);
        visTargets = unique(tVisTargets);
        visTargets = visTargets(2:end);        
        visHR = zeros(1,length(visTargets));
        for i = 1:length(visTargets)
            ind = tVisTargets == visTargets(i);
            visHR(i) = sum(ind & hit) ./ sum(ind & (hit | miss));
        end
        
        tAudTargets = celleqel2mat_padded(mworks.tSoundTargetAmplitude);
        audTargets = unique(tAudTargets);
        audTargets = audTargets(2:end);
        audHR = zeros(1,length(audTargets));
        for i = 1:length(audTargets)
            ind = tAudTargets == audTargets(i);
            audHR(i) = sum(ind & hit) ./ sum(ind & (hit | miss));
        end
        
        bxStruct(iexp).sn = mouseName;
        bxStruct(iexp).pctEarly = pctEarly;
        bxStruct(iexp).visHR = visHR;
        bxStruct(iexp).audHR = audHR;
        bxStruct(iexp).stimOnTime = mworks.stimOnTimeMs;
        bxStruct(iexp).stimOffTime = mworks.stimOffTimeMs;
        bxStruct(iexp).hit = hit;
        bxStruct(iexp).miss = miss;
        bxStruct(iexp).fa = fa;
        bxStruct(iexp).tVisTargets = tVisTargets;
        bxStruct(iexp).tAudTargets = tAudTargets;
        bxStruct(iexp).trLength = cell2mat(mworks.tCyclesOn);
        bxStruct(iexp).trLengthMs = cellfun(@(x,y) x-y,...
            mworks.tLeverReleaseTimeMs, mworks.tLeverPressTimeMs);
        
        if mworks.doShortCatchTrial
            if mworks.doOriDetect == 1
                visualBlock = 0;
                auditoryBlock = 1;
                bxStruct(iexp).invType = 'vis';
                bxStruct(iexp).tInvTargets = cell2mat_padded(...
                    mworks.tCatchGratingDirectionDeg);
            else
                visualBlock = 1;
                auditoryBlock = 0;
                bxStruct(iexp).invType = 'aud';
                bxStruct(iexp).tInvTargets = cell2mat_padded(...
                    mworks.tSoundCatchAmplitude);
            end
            bxStruct(iexp).invHit = strcmp(mworks.catchTrialOutcomeCell,'FA');
            bxStruct(iexp).invMiss = strcmp(mworks.catchTrialOutcomeCell,'CR');
            bxStruct(iexp).invTrLength = cell2mat(mworks.catchCyclesOn);
        else
            bxStruct(iexp).invType = nan;
            bxStruct(iexp).invHit = nan(1,nTrials);
            bxStruct(iexp).invMiss = nan(1,nTrials);
            bxStruct(iexp).invTrLength = nan(1,nTrials);
        end        
    end
end
