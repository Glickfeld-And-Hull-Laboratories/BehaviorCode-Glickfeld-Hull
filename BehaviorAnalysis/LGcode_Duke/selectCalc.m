function [s, b] = selectCalc(input,trials)
    
    if nargin == 2
        if ~isempty(trials)
            input = trialChopper(input,trials);
        end
    end
    if isfield(input,'doContrastDiscrim')
        if input.doContrastDiscrim
            rCon = celleqel2mat_padded(input.rightGratingContrast);
            lCon = celleqel2mat_padded(input.leftGratingContrast);
            if isempty(find(unique(rCon)==0))
                ratio_mat = chop(rCon./lCon,3);
                doRatio = 1;
            else
                ratio_mat = rCon-lCon;
                doRatio = 0;
            end
        elseif input.doSizeDiscrim
            rSize = celleqel2mat_padded(input.rightGratingDiameterDeg);
            lSize = celleqel2mat_padded(input.leftGratingDiameterDeg);
            ratio_mat = chop(rSize./lSize,3);
            doRatio = 1;
        elseif input.doOriDiscrim
            ratio_mat = -1.*celleqel2mat_padded(input.tGratingDirectionStart);
            doRatio = 0;
        end
    elseif isfield(input, 'rightGratingContrast')
        rCon = celleqel2mat_padded(input.rightGratingContrast);
        lCon = celleqel2mat_padded(input.leftGratingContrast);
        if ~find(unique(rCon)==0)
            ratio_mat = chop(rCon./lCon,3);
            doRatio = 1;
        else
            ratio_mat = rCon-lCon;
            doRatio = 0;
        end
        tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
    else
        tCon = celleqel2mat_padded(input.tGratingContrast);
        dCon = celleqel2mat_padded(input.dGratingContrast);
        tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
        rCon = zeros(size(tCon));
        lCon = zeros(size(tCon));
        rCon(find(tLeftTrial)) = dCon(find(tLeftTrial));
        rCon(find(~tLeftTrial)) = tCon(find(~tLeftTrial));
        lCon(find(~tLeftTrial)) = dCon(find(~tLeftTrial));
        lCon(find(tLeftTrial)) = tCon(find(tLeftTrial));
        if ~find(unique(rCon)==0)
            doRatio = 1;
            ratio_mat = chop(rCon./lCon,3);
        else
            ratio_mat = rCon-lCon;
            doRatio = 0;
        end
    end
    ratios = unique(ratio_mat);
    if isfield(input,'tRightResponse')
        tRight = celleqel2mat_padded(input.tRightResponse);
    else
        SIx = strcmp(input.trialOutcomeCell,'success');
        MIx = strcmp(input.trialOutcomeCell,'incorrect');
        tRight = zeros(size(rCon));
        tRight(intersect(find(~tLeftTrial),find(SIx))) = 1;
        tRight(intersect(find(tLeftTrial),find(MIx))) = 1;
    end
    if isfield(input,'doBlock2')
        if input.doBlock2
            b2ind = find(celleqel2mat_padded(input.tBlock2TrialNumber));
            ratio_mat(b2ind) = 0;
        end
    end
    if isfield(input,'isNoGo')
        indNogo = find(celleqel2mat_padded(input.isNoGo));
        ratio_mat(indNogo) = 0;
        tRight(indNogo) = 0;
    end
    Iind = find(strcmp(input.trialOutcomeCell,'ignore'));
    ratio_mat(Iind) = 0;
    totR = sum(ratio_mat==ratios(end));
    pctRonR = sum(ratio_mat==ratios(end) & tRight)./totR;
    totL = sum(ratio_mat==ratios(1));
    pctRonL = sum(ratio_mat==ratios(1) & tRight)./totL;
    s = pctRonR-pctRonL;
    if doRatio
        totR_all = sum(ratio_mat>1);
        totL_all = sum(ratio_mat<1 & ratio_mat>0);
        pctRonR_all = sum(ratio_mat>1 & tRight)./totR_all;
        pctRonL_all = sum(ratio_mat<1 & ratio_mat>0 & tRight)./totL_all;
    else
        totR_all = sum(ratio_mat>0);
        totL_all = sum(ratio_mat<0);
        pctRonR_all = sum(ratio_mat>0 & tRight)./totR_all;
        pctRonL_all = sum(ratio_mat<0 & tRight)./totL_all;
    end
    b = ((pctRonR_all+pctRonL_all)./2)-0.5;
end
