early_mat = zeros(nexp,1);
HR_max_mat = zeros(nexp,1);
bxExp = [];
for iexp = 1:nexp
    subnum = xd.Subject(iexp);
    expDate = char(xd.DateStr(iexp));
    
    mworks_dir = dir(fullfile(rc.pathStr, ['data-i' num2str(subnum) '-' expDate '-*']));
    runs = xd.ChooseMatFile(iexp);
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
    
%         nt = length(input.trialOutcomeCell);
%             try
%                 input_temp = [input_temp input];
%             catch
%                 inpNames1 = fieldnames(input_temp);
%                 inpNames2 = fieldnames(input);
%                 inpLong = gt(length(inpNames1),length(inpNames2));
%                 if inpLong == 1
%                     inpPlusInd = ismember(inpNames1,inpNames2);
%                     inpPlus = inpNames1(~inpPlusInd);
%                     for i = 1:length(inpPlus)
%                         input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
%                     end
%                 else
%                     inpPlusInd = ismember(inpNames2,inpNames1);
%                     inpPlus = inpNames2(~inpPlusInd);
%                     for i = 1:length(inpPlus)
%                         input_temp.(char(genvarname(inpPlus(i)))) = cell(1,nt);
%                     end
%                 end
%                 input_temp = [input_temp input];
%             end
%         end

    
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
    HR_max_mat(iexp) = pctCorr_maxOri;
    
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
    tOn = input_temp.stimOnTimeMs;
    tOff = input_temp.stimOffTimeMs;
    bxExp(iexp).tOn = tOn;
    bxExp(iexp).tOff = tOff;
    
    bxExp(iexp).sn = subnum;
end

save(fullfile(fnout,'bxExpMat_training'),'bxExp','early_mat','HR_max_mat')