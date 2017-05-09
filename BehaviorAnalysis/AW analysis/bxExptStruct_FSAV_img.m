early_mat = zeros(nexp,1);
HR_ori_mat = zeros(nexp,1);
HR_amp_mat = zeros(nexp,1);
bxExp = [];
for iexp = 1:nexp
    subnum = expt(iexp).SubNum;
    expDate = expt(iexp).date;
    nrun = expt(iexp).nrun;
    
    far_run = [];
    hr_run = [];
    input_temp = [];
    for irun = 1:nrun
        runTime = expt(iexp).time_mat(irun,:);
        mworks = ['data-i' subnum '-' expDate '-' runTime];
        load(fullfile(rc.pathStr,mworks))
        if irun == 1
            
            input_temp = input;   
            
        else
        nt = length(input.trialOutcomeCell);
            try
                input_temp = [input_temp input];
            catch
                inpNames1 = fieldnames(input_temp);
                inpNames2 = fieldnames(input);
                inpLong = gt(length(inpNames1),length(inpNames2));
                if inpLong == 1
                    inpPlusInd = ismember(inpNames1,inpNames2);
                    inpPlus = inpNames1(~inpPlusInd);
                    for i = 1:length(inpPlus)
                        input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
                    end
                else
                    inpPlusInd = ismember(inpNames2,inpNames1);
                    inpPlus = inpNames2(~inpPlusInd);
                    for i = 1:length(inpPlus)
                        input_temp.(char(genvarname(inpPlus(i)))) = cell(1,nt);
                    end
                end
                input_temp = [input_temp input];
            end
        end
        
        nt = length(input.trialOutcomeCell);
        far_temp = sum(strcmp(input.trialOutcomeCell,'failure'))./nt;
        tVisTar = round(cell2mat_padded(input.tGratingDirectionDeg));
        easyVis = max(unique(tVisTar));
        easyVis_ind = tVisTar == easyVis;
        easyVis_out = input.trialOutcomeCell(easyVis_ind);
        hr_temp = sum(strcmp(easyVis_out,'success'))./(sum(strcmp(easyVis_out,'success')) + sum(strcmp(easyVis_out,'ignore')));

        far_run(irun) = far_temp;
        hr_run(irun) = hr_temp;
    end
    bxExp(iexp).far_run = far_run;
    bxExp(iexp).hr_run = hr_run;
    
    input_temp = concatenateDataBlocks(input_temp);
    
    nt = length(input_temp.trialOutcomeCell);
    failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
    missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
    successIx = strcmp(input_temp.trialOutcomeCell, 'success');
    pctEarly = sum(failureIx,2)./length(failureIx);
    early_mat(iexp) = pctEarly;

    gratingDirectionDeg = cell2mat_padded(input_temp.tGratingDirectionDeg);
    if expt(iexp).catch
        bxExp(iexp).tInvVisTargets = cell2mat_padded(input_temp.tCatchGratingDirectionDeg);
    end
    oris = unique(gratingDirectionDeg);
    maxOriTrials = find(gratingDirectionDeg == max(oris,[],2));
    pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
    HR_ori_mat(iexp) = pctCorr_maxOri;
    
    bxExp(iexp).faIx = failureIx;
    bxExp(iexp).mIx = missedIx;
    bxExp(iexp).sIx = successIx;
    bxExp(iexp).tVisTargets = gratingDirectionDeg;
    
    if ~isfield(input_temp, 'tSoundTargetAmplitude')
        soundAmplitude = repmat(input_temp.soundTargetAmplitude,1,nt);
    else
        soundAmplitude = celleqel2mat_padded(input_temp.tSoundTargetAmplitude);
    end
    amps = unique(soundAmplitude);
    maxAmpTrials = find(soundAmplitude == max(amps,[],2));
    pctCorr_maxAmp = sum(successIx(maxAmpTrials),2)./(sum(successIx(maxAmpTrials),2)+sum(missedIx(maxAmpTrials),2));
    HR_amp_mat(iexp) = pctCorr_maxAmp;

    bxExp(iexp).tAudTargets = soundAmplitude;
    
    if ~isfield(input_temp, 'catchTrialOutcomeCell') | (sum(strcmp(input_temp.catchTrialOutcomeCell,'FA'))+sum(strcmp(input_temp.catchTrialOutcomeCell,'CR'))) == 0
        for trN = 1:length(input_temp.trialOutcomeCell)
            if input_temp.tShortCatchTrial{trN}
                if input_temp.tFalseAlarm{trN}
                    input_temp.catchTrialOutcomeCell{trN} = 'FA';
                end
                if isfield(input_temp, 'cCatchOn')
                    if isempty(input_temp.cCatchOn{trN})
                        input_temp.cCatchOn{trN} = NaN;
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
                    end
                    if (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})>input_temp.nFramesReact
                        input_temp.catchTrialOutcomeCell{trN} = 'CR';
                    end
                    if (input_temp.cLeverUp{trN}-input_temp.cCatchOn{trN})<input_temp.nFramesTooFast
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
                    end
                else
                    if isempty(input_temp.tCatchTimeMs{trN})
                        input_temp.cCatchOn{trN} = NaN;
                        input_temp.catchTrialOutcomeCell{trN} = 'failure';
                    end
                    if (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})>input_temp.reactTimeMs
                        input_temp.catchTrialOutcomeCell{trN} = 'CR';
                    end
                    if (input_temp.leverUpTimeMs{trN}-input_temp.tCatchTimeMs{trN})<input_temp.tooFastTimeMs
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
        
    bxExp(iexp).trLength = cell2mat(input_temp.tCyclesOn);
    if iscell(input_temp.nFramesOn);
        cOn = unique(cell2mat(input_temp.nFramesOn));
        cOff = unique(cell2mat(input_temp.nFramesOff));
    else
        cOn = input_temp.nFramesOn;
        cOff = input_temp.nFramesOff;
    end
    tOn = input_temp.stimOnTimeMs;
    tOff = input_temp.stimOffTimeMs;
    bxExp(iexp).cOn = cOn;
    bxExp(iexp).cOff = cOff;
    bxExp(iexp).tOn = tOn;
    bxExp(iexp).tOff = tOff;
    
    bxExp(iexp).sn = subnum;
end