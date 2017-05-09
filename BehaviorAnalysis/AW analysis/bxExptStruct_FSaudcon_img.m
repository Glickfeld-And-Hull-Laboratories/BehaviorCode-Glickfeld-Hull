early_mat = zeros(nexp,1);
HR_ori_mat = zeros(nexp,1);
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
    v_tr = cell2mat(input_temp.tBlock2TrialNumber) == 0;
    av_tr = cell2mat(input_temp.tBlock2TrialNumber) == 1;
    oris = unique(gratingDirectionDeg);
    maxOriTrials = find(gratingDirectionDeg == max(oris,[],2));
    pctCorr_maxOri = sum(successIx(maxOriTrials),2)./(sum(successIx(maxOriTrials),2)+sum(missedIx(maxOriTrials),2));
    HR_ori_mat(iexp) = pctCorr_maxOri;
    
    bxExp(iexp).allTargets = gratingDirectionDeg;
    bxExp(iexp).trType = cell2mat(input_temp.tBlock2TrialNumber);
    bxExp(iexp).vis_tr = v_tr;
    bxExp(iexp).audvis_tr = av_tr;
    bxExp(iexp).faIx = failureIx;
    bxExp(iexp).mIx = missedIx;
    bxExp(iexp).sIx = successIx;
    bxExp(iexp).tVisTargets = gratingDirectionDeg(v_tr);
    bxExp(iexp).tAVTargets = gratingDirectionDeg(av_tr);
        
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