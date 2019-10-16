close all
clear all
mouse = 'i433';
date_rng = '1907';
min_trials = 10;
s_min = .7;
behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
expt_mat = dir(fullfile(behav_path, ['data-' mouse '-' date_rng '*']));
t = 0;
date_mat = [];
clear temp
for i = 1:length(expt_mat) 
    load(fullfile(behav_path,expt_mat(i).name))
    if isfield(input,'leftGratingContrast') & input.doSizeDiscrim
        fprintf([expt_mat(i).name(11:16) '\n'])
        lCon = celleqel2mat_padded(input.leftGratingContrast);
        rCon = celleqel2mat_padded(input.rightGratingContrast);
        tCon = celleqel2mat_padded(input.tGratingContrast);
        dCon = celleqel2mat_padded(input.dGratingContrast);
        conRat_side = chop(rCon./lCon,2);
        conRat_targ = chop(tCon./dCon,2);
        if find(unique(conRat_side)>1)
            [s b] = selectCalc(input);
            if s > s_min & length(lCon)>100
                t = t+1;
                lSize = celleqel2mat_padded(input.leftGratingDiameterDeg);
                rSize = celleqel2mat_padded(input.rightGratingDiameterDeg);
                tSize = celleqel2mat_padded(input.tGratingDiameterDeg);
                dSize = celleqel2mat_padded(input.dGratingDiameterDeg);
                sizeRat_side = chop(rSize./lSize,2);
                sizeRat_targ = chop(tSize./dSize,2);
                rChoice = zeros(size(lCon));
                rChoice(intersect(find(sizeRat_side>1), find(strcmp(input.trialOutcomeCell,'success')))) = 1;
                rChoice(intersect(find(sizeRat_side<1), find(strcmp(input.trialOutcomeCell,'incorrect')))) = 1;
                pctCorr = strcmp(input.trialOutcomeCell,'success');
                size_side = unique(sizeRat_side);
                cons_side = unique(conRat_side);
                nSize_side = length(size_side(~isnan(size_side)));
                nCon_side = length(cons_side(~isnan(cons_side)));
                size_targ = unique(sizeRat_targ);
                cons_targ = unique(conRat_targ);
                nSize_targ = length(size_targ(~isnan(size_targ)));
                nCon_targ = length(cons_targ(~isnan(cons_targ)));
                pct_right_size_side = zeros(1,nSize_side);
                ci_size_side = zeros(2,nSize_side);
                for iSz = 1:nSize_side
                    ind = find(sizeRat_side == size_side(iSz));
                    [pct_right_size_side(1,iSz) ci_size_side(:,iSz)] = binofit(sum(rChoice(ind)),length(ind));
                end
                pct_right_con_side= zeros(1,nCon_side);
                ci_con_side = zeros(2,nCon_side);
                for iCon = 1:nCon_side
                    ind = find(conRat_side == cons_side(iCon));
                    [pct_right_con_side(1,iCon) ci_con_side(:,iCon)] = binofit(sum(rChoice(ind)),length(ind));
                end
                pct_right_size_targ = zeros(1,nSize_targ);
                ci_size_targ = zeros(2,nSize_targ);
                for iSz = 1:nSize_targ
                    ind = find(sizeRat_targ == size_targ(iSz));
                    [pct_right_size_targ(1,iSz) ci_size_targ(:,iSz)] = binofit(sum(pctCorr(ind)),length(ind));
                end
                pct_right_con_targ= zeros(1,nCon_targ);
                ci_con_targ = zeros(2,nCon_targ);
                for iCon = 1:nCon_targ
                    ind = find(conRat_targ == cons_targ(iCon));
                    [pct_right_con_targ(1,iCon) ci_con_targ(:,iCon)] = binofit(sum(pctCorr(ind)),length(ind));
                end
                figure;
                subplot(2,2,1)
                errorbar(size_side,pct_right_size_side, pct_right_size_side-ci_size_side(1,:),ci_size_side(2,:)-pct_right_size_side,'-o'); ax = gca; ax.XScale = 'log';
                xlabel('R/L size')
                ylabel('%Right choice')
                ylim([0 1])
                subplot(2,2,2)
                errorbar(cons_side,pct_right_con_side, pct_right_con_side-ci_con_side(1,:),ci_con_side(2,:)-pct_right_con_side,'-o'); ax = gca; ax.XScale = 'log';
                xlabel('R/L contrast')
                ylabel('%Right choice')
                ylim([0 1])
                subplot(2,2,3)
                errorbar(size_targ,pct_right_size_targ, pct_right_size_targ-ci_size_targ(1,:),ci_size_targ(2,:)-pct_right_size_targ,'-o'); ax = gca; ax.XScale = 'log';
                xlabel('T/D size')
                ylabel('%Correct')
                ylim([0 1])
                subplot(2,2,4)
                errorbar(cons_targ,pct_right_con_targ, pct_right_con_targ-ci_con_targ(1,:),ci_con_targ(2,:)-pct_right_con_targ,'-o'); ax = gca; ax.XScale = 'log';
                xlabel('T/D contrast')
                ylabel('%Correct')
                ylim([0 1])
                suptitle(expt_mat(i).name)
                temp(t).leftGratingContrast = input.leftGratingContrast;
                temp(t).rightGratingContrast = input.rightGratingContrast;
                temp(t).tGratingContrast = input.tGratingContrast;
                temp(t).dGratingContrast = input.dGratingContrast;
                temp(t).leftGratingDiameterDeg = input.leftGratingDiameterDeg;
                temp(t).rightGratingDiameterDeg = input.rightGratingDiameterDeg;
                temp(t).tGratingDiameterDeg = input.tGratingDiameterDeg;
                temp(t).dGratingDiameterDeg = input.dGratingDiameterDeg;
                temp(t).trialOutcomeCell = input.trialOutcomeCell;
                date_mat = [date_mat; expt_mat(i).name(11:16)];
            end
        end
    end
end
if size(date_mat,1)>0
    all_trials = concatenateDataBlocks(temp);
    lCon = celleqel2mat_padded(all_trials.leftGratingContrast);
    rCon = celleqel2mat_padded(all_trials.rightGratingContrast);
    tCon = celleqel2mat_padded(all_trials.tGratingContrast);
    dCon = celleqel2mat_padded(all_trials.dGratingContrast);
    conRat_side = chop(rCon./lCon,2);
    conRat_targ = chop(tCon./dCon,2);
    lSize = celleqel2mat_padded(all_trials.leftGratingDiameterDeg);
    rSize = celleqel2mat_padded(all_trials.rightGratingDiameterDeg);
    tSize = celleqel2mat_padded(all_trials.tGratingDiameterDeg);
    dSize = celleqel2mat_padded(all_trials.dGratingDiameterDeg);
    sizeRat_side = chop(rSize./lSize,2);
    sizeRat_targ = chop(tSize./dSize,2);
    rChoice = zeros(size(lCon));
    rChoice(intersect(find(sizeRat_side>1), find(strcmp(all_trials.trialOutcomeCell,'success')))) = 1;
    rChoice(intersect(find(sizeRat_side<1), find(strcmp(all_trials.trialOutcomeCell,'incorrect')))) = 1;
    pctCorr = strcmp(all_trials.trialOutcomeCell,'success');
    size_side = unique(sizeRat_side);
    cons_side = unique(conRat_side);
    nSize_side = length(size_side);
    nCon_side = length(cons_side);
    size_targ = unique(sizeRat_targ);
    cons_targ = unique(conRat_targ);
    nSize_targ = length(size_targ);
    nCon_targ = length(cons_targ);
    pct_right_size_side = nan(1,nSize_side);
    ci_size_side = nan(2,nSize_side);
    for iSz = 1:nSize_side
        ind = find(sizeRat_side == size_side(iSz));
        if length(ind)>min_trials
            [pct_right_size_side(1,iSz) ci_size_side(:,iSz)] = binofit(sum(rChoice(ind)),length(ind));
        end
    end
    pct_right_con_side= nan(1,nCon_side);
    ci_con_side = nan(2,nCon_side);
    for iCon = 1:nCon_side
        ind = find(conRat_side == cons_side(iCon));
        if length(ind)>min_trials
            [pct_right_con_side(1,iCon) ci_con_side(:,iCon)] = binofit(sum(rChoice(ind)),length(ind));
        end
    end
    pct_right_size_targ = nan(1,nSize_targ);
    ci_size_targ = nan(2,nSize_targ);
    for iSz = 1:nSize_targ
        ind = find(sizeRat_targ == size_targ(iSz));
        if length(ind)>min_trials
        [pct_right_size_targ(1,iSz) ci_size_targ(:,iSz)] = binofit(sum(pctCorr(ind)),length(ind));
        end
    end
    pct_right_con_targ= nan(1,nCon_targ);
    ci_con_targ = nan(2,nCon_targ);
    for iCon = 1:nCon_targ
        ind = find(conRat_targ == cons_targ(iCon));
        if length(ind)>min_trials
            [pct_right_con_targ(1,iCon) ci_con_targ(:,iCon)] = binofit(sum(pctCorr(ind)),length(ind));
        end
    end
    figure;
    subplot(2,2,1)
    errorbar(size_side,pct_right_size_side, pct_right_size_side-ci_size_side(1,:),ci_size_side(2,:)-pct_right_size_side,'-o'); ax = gca; ax.XScale = 'log';
    xlabel('R/L size')
    ylabel('%Right choice')
    ylim([0 1])
    subplot(2,2,2)
    errorbar(cons_side,pct_right_con_side, pct_right_con_side-ci_con_side(1,:),ci_con_side(2,:)-pct_right_con_side,'-o'); ax = gca; ax.XScale = 'log';
    xlabel('R/L contrast')
    ylabel('%Right choice')
    ylim([0 1])
    subplot(2,2,3)
    errorbar(size_targ,pct_right_size_targ, pct_right_size_targ-ci_size_targ(1,:),ci_size_targ(2,:)-pct_right_size_targ,'-o'); ax = gca; ax.XScale = 'log';
    xlabel('T/D size')
    ylabel('%Correct')
    ylim([0 1])
    subplot(2,2,4)
    errorbar(cons_targ,pct_right_con_targ, pct_right_con_targ-ci_con_targ(1,:),ci_con_targ(2,:)-pct_right_con_targ,'-o'); ax = gca; ax.XScale = 'log';
    xlabel('T/D contrast')
    ylabel('%Correct')
    ylim([0 1])
    suptitle([mouse ': ' date_mat(1,:) '-' date_mat(end,:) '- ' num2str(size(date_mat,1)) ' sessions; ' num2str(size(conRat_targ,2)) ' trials'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' mouse '_Summary_SizeVCon.pdf'],'-dpdf','-bestfit'); 
    save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' mouse '_Summary_SizeVCon.mat'],'all_trials')
end