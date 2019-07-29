function [s, b] = selectCalc(input,trials)
    
    if nargin == 2
        if ~isempty(trials)
            input = trialChopper(input,trials);
        end
    end
    
    if input.doContrastDiscrim
        rCon = celleqel2mat_padded(input.rightGratingContrast);
        lCon = celleqel2mat_padded(input.leftGratingContrast);
        ratio_mat = chop(rCon./lCon,3);
    elseif input.doSizeDiscrim
        rSize = celleqel2mat_padded(input.rightGratingDiameterDeg);
        lSize = celleqel2mat_padded(input.leftGratingDiameterDeg);
        ratio_mat = chop(rSize./lSize,3);
    elseif input.doOriDiscrim
        ratio_mat = celleqel2mat_padded(input.tGratingDirectionStart);
    end
    ratios = unique(ratio_mat);
    tRight = celleqel2mat_padded(input.tRightResponse);
    if input.doBlock2
        b2ind = find(celleqel2mat_padded(input.tBlock2TrialNumber));
        ratio_mat(b2ind) = 0;
    end
    Iind = find(strcmp(input.trialOutcomeCell,'ignore'));
    ratio_mat(Iind) = 0;
    totR = sum(ratio_mat==ratios(end));
    pctRonR = sum(ratio_mat==ratios(end) & tRight)./totR;
    totL = sum(ratio_mat==ratios(1));
    pctRonL = sum(ratio_mat==ratios(1) & tRight)./totL;
    s = pctRonR-pctRonL;
    totR_all = sum(ratio_mat>1);
    totL_all = sum(ratio_mat<1 & ratio_mat>0);
    pctRonR_all = sum(ratio_mat>1 & tRight)./totR_all;
    pctRonL_all = sum(ratio_mat<1 & ratio_mat>0 & tRight)./totL_all;
    b = ((pctRonR_all+pctRonL_all)./2)-0.5;
end
