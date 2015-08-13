function mouse = createExptStructAttn2015;
    rc = behavConstsAttn;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAttn;
    
    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).posN = [];
        mouse(imouse).taskN = [];
        mouse(imouse).fitN = [];
        for itask = 1:2
            mouse(imouse).task(itask).ind = [];
            for ipos = 1:2  
                mouse(imouse).task(itask).pos(ipos).ind = [];
                for iblock = 1:2
                    mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect = [];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed = [];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).testCon = [];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).b2RewardPct = [];
                    mouse(imouse).task(itask).block(iblock).nCorrect = [];
                    mouse(imouse).task(itask).block(iblock).nMissed = [];
                    mouse(imouse).task(itask).block(iblock).testCon = [];
                    mouse(imouse).task(itask).block(iblock).b2RewardPct = [];
                    for ihalf = 1:2
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect = [];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed = [];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).testCon = [];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).b2RewardPct = [];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect = [];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed = [];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).testCon = [];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).b2RewardPct = [];
                    end
                end
            end
        end
        for iexp = 1:length(mouse(imouse).ind)
            for iblock = 1:2
                mouse(imouse).expt(iexp).block(iblock).nCorrect = 0;
                mouse(imouse).expt(iexp).block(iblock).nMissed = 0;
                for ihalf = 1:2
                    mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect = 0;
                    mouse(imouse).expt(iexp).block(iblock).half(ihalf).nMissed = 0;
                end
            end
            ind = mouse(imouse).ind(iexp);
            fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(ind)) '-' cell2mat(xd.DateStr(ind)) '.mat']);
            load(fn)
            
            trials = str2num(cell2mat(xd.TrialRangeToUse(ind)));
            if isnan(trials)
                trials = input.trialSinceReset;
                range = 1:trials;
            else
                range = trials;
            end
            midpoint = ceil((range(end)-range(1))/2);
            firsthalftrials = range(1):midpoint;
            secondhalftrials = midpoint+1:range(end);
            mouse(imouse).expt(iexp).block(1).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==0), range);
            mouse(imouse).expt(iexp).block(2).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==1), range);
            mouse(imouse).expt(iexp).block(1).intensities = unique(chop(unique(cell2mat(input.gratingContrast)),2));                
            mouse(imouse).expt(iexp).block(2).intensities = unique(chop(unique(cell2mat(input.gratingContrast(mouse(imouse).expt(iexp).block(2).trial_ind))),2));                    
            if isfield(input, 'doBlock2SeparateReward')
                if input.doBlock2SeparateReward
                    mouse(imouse).expt(iexp).doB2reward = 1;
                    mouse(imouse).expt(iexp).reward = [(input.maxRewardUs+input.maxRewardUs)./2 (input.block2MaxRewardUs+input.block2MaxRewardUs)./2];   
                    mouse(imouse).expt(iexp).b2RewardPct = ((input.block2MaxRewardUs+input.block2MaxRewardUs)./2)./((input.maxRewardUs+input.maxRewardUs)./2);
                else
                    mouse(imouse).expt(iexp).doB2reward = 0;
                    mouse(imouse).expt(iexp).reward = [(input.maxRewardUs+input.maxRewardUs)./2 (input.maxRewardUs+input.maxRewardUs)./2];
                    mouse(imouse).expt(iexp).b2RewardPct = 1;
                end
            else
                mouse(imouse).expt(iexp).doB2reward = 0;
                mouse(imouse).expt(iexp).reward = [(input.maxRewardUs+input.maxRewardUs)./2 (input.maxRewardUs+input.maxRewardUs)./2];   
                mouse(imouse).expt(iexp).b2RewardPct = 1;
            end
            if isfield(input, 'doBlock2SeparateOdds')
                if input.doBlock2SeparateOdds
                    mouse(imouse).expt(iexp).doB2Odds = 1;
                    mouse(imouse).expt(iexp).b2OddsPct = sum(input.block2TrPer80V,2)./80;
                    mouse(imouse).expt(iexp).odds = input.block2TrPer80V;   
                else
                    mouse(imouse).expt(iexp).odds = input.trPer80V;
                    mouse(imouse).expt(iexp).doB2Odds = 0;
                    mouse(imouse).expt(iexp).b2OddsPct = .5;
                end
            else
                mouse(imouse).expt(iexp).odds = input.trPer80V;
                mouse(imouse).expt(iexp).doB2Odds = 0;
                mouse(imouse).expt(iexp).b2OddsPct = .5;
            end
             if input.gratingAzimuthDeg > 0
                ipos = 1;
            else
                ipos = 2;
            end
            if mouse(imouse).expt(iexp).doB2reward == 0
                itask = 1;
            else
                itask = 2;
            end
            mouse(imouse).expt(iexp).pos = ipos;
            mouse(imouse).expt(iexp).task = itask;
            mouse(imouse).posN = [mouse(imouse).posN mouse(imouse).expt(iexp).pos];
            mouse(imouse).taskN = [mouse(imouse).taskN mouse(imouse).expt(iexp).task];
            fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(ind)) '-max.mat']);
            load(fn);
            mouse(imouse).expt(iexp).block(1).thresh = fitS(1).thresh;
            mouse(imouse).expt(iexp).block(1).ci95 = fitS(1).bootStats.ci95;
            ci95_ratio_1 = fitS(1).bootStats.ci95(2)./fitS(1).bootStats.ci95(1);
            if ci95_ratio_1 < 2
                mouse(imouse).expt(iexp).well_fit = 1;
                mouse(imouse).expt(iexp).intensityX = fitS(1).intensityX;
                mouse(imouse).expt(iexp).fractCorrY = fitS(1).fractCorrY;
                mouse(imouse).expt(iexp).coefEsts = fitS(1).coefEsts;
                mouse(imouse).task(itask).pos(ipos).ind = [mouse(imouse).task(itask).pos(ipos).ind iexp];
                mouse(imouse).task(itask).ind = [mouse(imouse).task(itask).ind iexp];
                mouse(imouse).fitN = [mouse(imouse).fitN 1];
            else
                mouse(imouse).expt(iexp).block(1).well_fit = 0;
                mouse(imouse).fitN = [mouse(imouse).fitN 0];
            end
            for iblock = 1:2
                for icon = 1:length(mouse(imouse).expt(iexp).block(2).intensities)
                    con_ind = find(chop(cell2mat(input.gratingContrast),2)==mouse(imouse).expt(iexp).block(2).intensities(icon));
                    mouse(imouse).expt(iexp).block(iblock).nTrials(icon) = length(intersect(mouse(imouse).expt(iexp).block(iblock).trial_ind, con_ind));
                    mouse(imouse).expt(iexp).block(iblock).nCorrect(icon) = length(intersect(mouse(imouse).expt(iexp).block(iblock).trial_ind, intersect(con_ind, intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                    mouse(imouse).expt(iexp).block(iblock).nMissed(icon) = length(intersect(mouse(imouse).expt(iexp).block(iblock).trial_ind, intersect(con_ind, find(cell2mat(input.reactTimesMs)>550))));
                    [x y] = binofit(mouse(imouse).expt(iexp).block(iblock).nCorrect(icon),mouse(imouse).expt(iexp).block(iblock).nCorrect(icon)+mouse(imouse).expt(iexp).block(iblock).nMissed(icon));
                    mouse(imouse).expt(iexp).block(iblock).binofit(icon).pctCorrect = x;
                    mouse(imouse).expt(iexp).block(iblock).binofit(icon).ci95 = y;
                    for ihalf = 1:2
                        if ihalf  == 1
                            half_trials = intersect(mouse(imouse).expt(iexp).block(iblock).trial_ind,firsthalftrials);
                        elseif ihalf  == 2
                            half_trials = intersect(mouse(imouse).expt(iexp).block(iblock).trial_ind,secondhalftrials);
                        end
                        mouse(imouse).expt(iexp).block(iblock).half(ihalf).nTrials(icon) = length(intersect(half_trials, con_ind));
                        mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect(icon) = length(intersect(half_trials, intersect(con_ind, intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                        mouse(imouse).expt(iexp).block(iblock).half(ihalf).nMissed(icon) = length(intersect(half_trials, intersect(con_ind, find(cell2mat(input.reactTimesMs)>550))));
                        [x y] = binofit(mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect(icon),mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect(icon)+mouse(imouse).expt(iexp).block(iblock).half(ihalf).nMissed(icon));
                        mouse(imouse).expt(iexp).block(iblock).binofit(icon).half(ihalf).pctCorrect = x;
                        mouse(imouse).expt(iexp).block(iblock).binofit(icon).half(ihalf).ci95 = y;
                    end
                end
                if ci95_ratio_1 < 2
                    mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect = [mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect mouse(imouse).expt(iexp).block(iblock).nCorrect(1)];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed = [mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed mouse(imouse).expt(iexp).block(iblock).nMissed(1)];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).testCon = [mouse(imouse).task(itask).pos(ipos).block(iblock).testCon mouse(imouse).expt(iexp).block(2).intensities(1)];
                    mouse(imouse).task(itask).pos(ipos).block(iblock).b2RewardPct = [mouse(imouse).task(itask).pos(ipos).block(iblock).b2RewardPct mouse(imouse).expt(iexp).b2RewardPct];
                    mouse(imouse).task(itask).block(iblock).nCorrect = [mouse(imouse).task(itask).block(iblock).nCorrect mouse(imouse).expt(iexp).block(iblock).nCorrect(1)];
                    mouse(imouse).task(itask).block(iblock).nMissed = [mouse(imouse).task(itask).block(iblock).nMissed mouse(imouse).expt(iexp).block(iblock).nMissed(1)];
                    mouse(imouse).task(itask).block(iblock).testCon = [mouse(imouse).task(itask).block(iblock).testCon mouse(imouse).expt(iexp).block(2).intensities(1)];
                    mouse(imouse).task(itask).block(iblock).b2RewardPct = [mouse(imouse).task(itask).block(iblock).b2RewardPct mouse(imouse).expt(iexp).b2RewardPct];
                    for ihalf = 1:2
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect = [mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect(1)];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed = [mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed mouse(imouse).expt(iexp).block(iblock).half(ihalf).nMissed(1)];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).testCon = [mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).testCon mouse(imouse).expt(iexp).block(2).intensities(1)];
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).b2RewardPct = [mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).b2RewardPct mouse(imouse).expt(iexp).b2RewardPct];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect = [mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect mouse(imouse).expt(iexp).block(iblock).half(ihalf).nCorrect(1)];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed = [mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed mouse(imouse).expt(iexp).block(iblock).half(ihalf).nMissed(1)];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).testCon = [mouse(imouse).task(itask).block(iblock).half(ihalf).testCon mouse(imouse).expt(iexp).block(2).intensities(1)];
                        mouse(imouse).task(itask).block(iblock).half(ihalf).b2RewardPct = [mouse(imouse).task(itask).block(iblock).half(ihalf).b2RewardPct mouse(imouse).expt(iexp).b2RewardPct];
                    end
                end
            end
        end
        for itask = 1:2
            for ipos = 1:2
                for iblock = 1:2
                    [x y] = binofit(sum(mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect,2),sum(mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect,2) + sum(mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed,2));
                    mouse(imouse).task(itask).pos(ipos).block(iblock).pctCorrect = x;
                    mouse(imouse).task(itask).pos(ipos).block(iblock).ci95 = y;
                    [x y] = binofit(sum(mouse(imouse).task(itask).block(iblock).nCorrect,2),sum(mouse(imouse).task(itask).block(iblock).nCorrect,2) + sum(mouse(imouse).task(itask).block(iblock).nMissed,2));
                    mouse(imouse).task(itask).block(iblock).pctCorrect = x;
                    mouse(imouse).task(itask).block(iblock).ci95 = y;
                    for ihalf = 1:2
                        [x y] = binofit(sum(mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect,2), sum(mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect,2) + sum(mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed,2));
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).pctCorrect = x;
                        mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).ci95 = y;
                        [x y] = binofit(sum(mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect,2), sum(mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect,2) + sum(mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed,2));
                        mouse(imouse).task(itask).block(iblock).half(ihalf).pctCorrect = x;
                        mouse(imouse).task(itask).block(iblock).half(ihalf).ci95 = y;
                    end
                end
            end
        end
        mouse(imouse).pos_diff = diff(mouse(imouse).posN);
        mouse(imouse).modelFun = fitS(1).modelFun;
    end
    save([rc.fitOutputSummary '\' date '-attn_mouse.mat'], 'mouse')
end

    
    