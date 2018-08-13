stimOnTime100 = zeros(nexp,1);
early_mat = zeros(nexp,1);
HR_ori_mat = zeros(nexp,1);
HR_amp_mat = zeros(nexp,1);
visHitRate = cell(1,nexp);
audHitRate = cell(1,nexp);
bxExp = [];
for iexp = 1:nexp
    fprintf([num2str(iexp) ' '])
    subnum = xd.Subject(iexp);
    expDate = char(xd.DateStr(iexp));    
    
    mworks_dir = dir(fullfile(rc.pathStr, ['data-i' num2str(subnum) '-' expDate '-*']));
    runs = xd.ChooseMatFile(iexp);
    if iscell(runs)
        runs = cell2mat(runs);
        if ischar(runs)
            if strcmp(runs,'NaN')
                runs = nan;
            else
                runs = eval(runs);
            end
        end
    end
    if ~isnan(runs)
        mworks_dir = mworks_dir(runs);
        
        input_temp = [];
        for irun = 1:length(runs)
            load(fullfile(rc.pathStr,mworks_dir(irun).name))
            input_temp = [input_temp, input];
        end
        input_temp = concatenateDataBlocks(input_temp);
    else
        load(fullfile(rc.pathStr,mworks_dir(1).name))
        input_temp = input;
    end
    rng = str2num(xd.TrialRangeToUse{(iexp)});
    if ~isnan(rng(1))
        if rng(end) > input_temp.trialSinceReset
            input_temp = trialChopper(input_temp, [rng(1) input_temp.trialSinceReset]);
        else
            input_temp = trialChopper(input_temp, [rng(1) rng(end)]);
        end
    end
    
    stimOnTime100(iexp) = input_temp.stimOnTimeMs == 100;
    
    nt = length(input_temp.trialOutcomeCell);
    failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
    missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
    successIx = strcmp(input_temp.trialOutcomeCell, 'success');
    pctEarly = sum(failureIx,2)./length(failureIx);
    early_mat(iexp) = pctEarly;

    if input.doShortCatchTrial
        if input.doOriDetect == 1
            vis = 0;
            aud = 1;
            bxExp(iexp).invType = 'vis';
            bxExp(iexp).tInvTargets = cell2mat_padded(input_temp.tCatchGratingDirectionDeg);
        else
            vis = 1;
            aud = 0;
            bxExp(iexp).invType = 'aud';
            bxExp(iexp).tInvTargets = celleqel2mat_padded(input_temp.tSoundCatchAmplitude);
        end
    else
        bxExp(iexp).invType = NaN;
    end
    
    gratingDirectionDeg = cell2mat_padded(input_temp.tGratingDirectionDeg);
    oris = unique(gratingDirectionDeg);
    maxOriTrials = find(gratingDirectionDeg == max(oris,[],2));
    pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
    HR_ori_mat(iexp) = pctCorr_maxOri;
    
    
    if ~isfield(input_temp, 'tSoundTargetAmplitude')
        soundAmplitude = repmat(input_temp.soundTargetAmplitude,1,nt);
    else
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);
    end
    amps = unique(soundAmplitude);
    maxAmpTrials = find(soundAmplitude == max(amps,[],2));
    pctCorr_maxAmp = sum(successIx(maxAmpTrials),2)./(sum(successIx(maxAmpTrials),2)+sum(missedIx(maxAmpTrials),2));
    HR_amp_mat(iexp) = pctCorr_maxAmp;

    bxExp(iexp).faIx = failureIx;
    bxExp(iexp).mIx = missedIx;
    bxExp(iexp).sIx = successIx;
    bxExp(iexp).tVisTargets = gratingDirectionDeg;
    bxExp(iexp).tAudTargets = soundAmplitude;
    
    visHitRate{iexp} = zeros(1,length(oris)-1);
    for i = 1:length(oris)-1
        ind = gratingDirectionDeg == oris(i+1);
        visHitRate{iexp}(i) = sum(ind & successIx)./sum(ind & (successIx | missedIx));
    end
    audHitRate{iexp} = zeros(1,length(amps)-1);
    for i = 1:length(amps)-1
        ind = soundAmplitude == amps(i+1);
        audHitRate{iexp}(i) = sum(ind & successIx)./sum(ind & (successIx | missedIx));
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
    
    bxExp(iexp).invHitIx = invHit;
    bxExp(iexp).invMissIx = invMiss;
        
    
    leverUpTime = cell2mat(input_temp.leverUpTimeMs);
    leverDownTime = cell2mat(input_temp.leverDownTimeMs);
    trialTimeMs = double(leverUpTime - leverDownTime);

    tOn = double(input_temp.stimOnTimeMs);
    tOff = double(input_temp.stimOffTimeMs);
    
%     tInvCyc = cell2mat(input_temp.catchCyclesOn);
    catchTimeMs = celleqel2mat_padded(input_temp.tCatchTimeMs) - double(leverDownTime);
    catchCycCalc = round(catchTimeMs./double(tOn+tOff))+1;
    
    nCyc = cell2mat(input_temp.nCyclesOn);
    
    invTargetOnTime = double((tOn+tOff).*catchCycCalc);
    valTargetOnTime = double((tOn+tOff).*nCyc);

    valReactTimeCalc = trialTimeMs - valTargetOnTime;
    invReactTimeCalc = trialTimeMs - invTargetOnTime;
    
    bxExp(iexp).trLength = trialTimeMs;%cell2mat(input_temp.tCyclesOn);
    bxExp(iexp).nCyc = nCyc;
    bxExp(iexp).invTrLength = catchCycCalc;
    bxExp(iexp).valReact = valReactTimeCalc;
    bxExp(iexp).invReact = invReactTimeCalc;
    bxExp(iexp).tOn = tOn;
    bxExp(iexp).tOff = tOff;
    
    bxExp(iexp).sn = subnum;
end

save(fullfile(fnout,'bxExpMat_training'),'bxExp',...
'early_mat','HR_ori_mat','HR_amp_mat','visHitRate','audHitRate','stimOnTime100')