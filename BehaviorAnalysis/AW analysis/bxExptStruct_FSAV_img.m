nexp_img = size(imagingDatasets,2);
stimOnTime100_img = zeros(nexp_img,1);
early_mat_img = zeros(nexp_img,1);
HR_ori_mat_img = zeros(nexp_img,1);
HR_amp_mat_img = zeros(nexp_img,1);
visHitRate_img = cell(1,nexp_img);
audHitRate_img = cell(1,nexp_img);
bxImgExp = [];
for iexp = 1:nexp_img
    fprintf([num2str(iexp) ' '])
    subnum = imagingDatasets(iexp).SubNum;
    expDate = imagingDatasets(iexp).date;    
    
    mworks_dir = fullfile(rc.pathStr, ['data-i' num2str(subnum) '-' expDate '-']);
    runs = num2cell(imagingDatasets(iexp).time_mat,2);
    if length(runs) > 1        
        input_temp = [];
        for irun = 1:length(runs)
            load([mworks_dir runs{irun}])
            try
                input_temp = [input_temp, input];
            catch
                inpFields1 = fieldnames(input_temp);
                inpFields2 = fieldnames(input);
                inpFieldsPlus = unique(cat(1,setdiff(inpFields1,inpFields2),...
                    setdiff(inpFields2,inpFields1)));
                for i = 1:length(inpFieldsPlus)
                    if ~isfield(input_temp,inpFieldsPlus{i})
                        input_temp.(genvarname(inpFieldsPlus{i})) = ...
                            cell(1,input_temp.trialSinceReset);
                    elseif ~isfield(input,inpFieldsPlus{i})
                        input.(genvarname(inpFieldsPlus{i})) = ...
                            cell(1,input.trialSinceReset);
                    end
                end
                input_temp = [input_temp, input];
            end
        end
        input_temp = concatenateDataBlocks(input_temp);
    else
        load([mworks_dir runs{1}])
        input_temp = input;
    end
    stimOnTime100_img(iexp) = input_temp.stimOnTimeMs == 100;
    nt = length(input_temp.trialOutcomeCell);
    failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
    missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
    successIx = strcmp(input_temp.trialOutcomeCell, 'success');
    gratingDirectionDeg = cell2mat_padded(input_temp.tGratingDirectionDeg);
    if ~isfield(input_temp, 'tSoundTargetAmplitude')
        soundAmplitude = repmat(input_temp.soundTargetAmplitude,1,nt);
    else
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);
    end
    if input.doShortCatchTrial
        if input.doOriDetect == 1
            vis = 0;
            aud = 1;
            bxImgExp(iexp).invType = 'vis';
            bxImgExp(iexp).tInvTargets = cell2mat_padded(input_temp.tCatchGratingDirectionDeg);
        else
            vis = 1;
            aud = 0;
            bxImgExp(iexp).invType = 'aud';
            bxImgExp(iexp).tInvTargets = cell2mat_padded(input_temp.tSoundCatchAmplitude);
        end
    else
        bxImgExp(iexp).invType = NaN;
    end
    if ~isfield(input_temp, 'catchTrialOutcomeCell') | (sum(strcmp(input_temp.catchTrialOutcomeCell,'FA'))+sum(strcmp(input_temp.catchTrialOutcomeCell,'CR'))) == 0
        input_temp.catchTrialOutcomeCell = cell(1,length(gratingDirectionDeg));
        for trN = 1:length(input_temp.trialOutcomeCell)
            if input_temp.tShortCatchTrial{trN}
                if input_temp.tFalseAlarm{trN}
                    input_temp.catchTrialOutcomeCell{trN} = 'FA';
%                 end
                elseif isfield(input_temp, 'cCatchOn')
                    if isempty(input_temp.cCatchOn{trN})
                        input_temp.cCatchOn{trN} = NaN;
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
%                     end
                    elseif (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})>input_temp.nFramesReact
                        input_temp.catchTrialOutcomeCell{trN} = 'CR';
%                     end
                    elseif (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})<input_temp.nFramesTooFast
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
                    end
                else
                    if isempty(input_temp.tCatchTimeMs{trN})
                        input_temp.cCatchOn{trN} = NaN;
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
%                     end
                    elseif (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})>input_temp.reactTimeMs
                        input_temp.catchTrialOutcomeCell{trN} = 'CR';
%                     end
                    elseif (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})<input_temp.tooFastTimeMs
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
                    end
                end
            else
                input_temp.catchTrialOutcomeCell{trN} = 'NaN';
            end
        end
    end
    invHit = strcmp(input_temp.catchTrialOutcomeCell,'FA');
    invMiss = strcmp(input_temp.catchTrialOutcomeCell,'CR');
    trLength = cell2mat(input_temp.tCyclesOn);
%     invTrLength = cell2mat(input_temp.catchCyclesOn); 
    nCyc = cell2mat(input_temp.nCyclesOn);
    
    leverUpTime = cell2mat(input_temp.cLeverUp);
    leverDownTime = cell2mat(input_temp.cLeverDown);
    trialTimeFrames = leverUpTime - leverDownTime;
    trialTimeMs = double(trialTimeFrames .* (1000./input_temp.frameRateHz));
    
    catchTimeFrames = celleqel2mat_padded(input_temp.cCatchOn) - double(leverDownTime);
    catchTimeMs = catchTimeFrames .* (1000./input_temp.frameRateHz);
    if ~isnan(imagingDatasets(iexp).trial_range)
        trial_range = imagingDatasets(iexp).trial_range;
        nt = length(trial_range);
        failureIx = failureIx(trial_range);
        
        missedIx = missedIx(trial_range);
        successIx = successIx(trial_range);
        gratingDirectionDeg = gratingDirectionDeg(trial_range);
        soundAmplitude = soundAmplitude(trial_range);
        if input.doShortCatchTrial
            bxImgExp(iexp).tInvTargets = bxImgExp(iexp).tInvTargets(trial_range);
        end
        invHit = invHit(trial_range);
        invMiss = invMiss(trial_range);
        trLength = trLength(trial_range);  
%         invTrLength = invTrLength(trial_range);
        nCyc = nCyc(trial_range);
        trialTimeMs = trialTimeMs(trial_range);
        catchTimeMs = catchTimeMs(trial_range);
    end
    pctEarly = sum(failureIx,2)./length(failureIx);
    early_mat_img(iexp) = pctEarly;
    
    
    oris = unique(gratingDirectionDeg);
    maxOriTrials = find(gratingDirectionDeg == max(oris,[],2));
    pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
    HR_ori_mat_img(iexp) = pctCorr_maxOri;
    
    visHitRate_img{iexp} = zeros(1,length(oris)-1);
    for i = 1:length(oris)-1
        ind = gratingDirectionDeg == oris(i+1);
        visHitRate_img{iexp}(i) = sum(ind & successIx)./sum(ind & (successIx | missedIx));
    end   
    
    amps = unique(soundAmplitude);
    maxAmpTrials = find(soundAmplitude == max(amps,[],2));
    pctCorr_maxAmp = sum(successIx(maxAmpTrials),2)./(sum(successIx(maxAmpTrials),2)+sum(missedIx(maxAmpTrials),2));
    HR_amp_mat_img(iexp) = pctCorr_maxAmp;
    
    audHitRate_img{iexp} = zeros(1,length(amps)-1);
    for i = 1:length(amps)-1
        ind = soundAmplitude == amps(i+1);
        audHitRate_img{iexp}(i) = sum(ind & successIx)./sum(ind & (successIx | missedIx));
    end

    bxImgExp(iexp).faIx = failureIx;
    bxImgExp(iexp).mIx = missedIx;
    bxImgExp(iexp).sIx = successIx;
    bxImgExp(iexp).tVisTargets = gratingDirectionDeg;
    bxImgExp(iexp).tAudTargets = soundAmplitude;
    
    bxImgExp(iexp).invHitIx = invHit;
    bxImgExp(iexp).invMissIx = invMiss;
        
    
    tOn = input_temp.stimOnTimeMs;
    tOff = input_temp.stimOffTimeMs;
    
    catchCycCalc = round(catchTimeMs/double(tOn+tOff));
    
    invTargetOnTime = double((tOn+tOff).*catchCycCalc);
    valTargetOnTime = double((tOn+tOff).*nCyc);

    valReactTimeCalc = double(trialTimeMs - valTargetOnTime);
    invReactTimeCalc = double(trialTimeMs - invTargetOnTime);
    
    bxImgExp(iexp).trLength = nCyc;
    bxImgExp(iexp).invTrLength = catchCycCalc;
    bxImgExp(iexp).valReact = valReactTimeCalc;
    bxImgExp(iexp).invReact = invReactTimeCalc;
    bxImgExp(iexp).tOn = tOn;
    bxImgExp(iexp).tOff = tOff;
    
    bxImgExp(iexp).sn = subnum;
end

save(fullfile(fnout,'bxExpMat_imaging'),...
    'bxImgExp','early_mat_img','HR_ori_mat_img','HR_amp_mat_img',...
    'visHitRate_img','audHitRate_img','stimOnTime100_img')