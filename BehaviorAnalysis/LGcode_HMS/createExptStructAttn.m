function mouse = createExptStructAttn(itask)
%takes all data from xls, sorts according to task type
    rc = behavConstsAttn;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAttn;

    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).task(itask).ind = intersect(find(xd.Task == itask), mouse(imouse).ind);
        for ipos = 1:length(pv.pos_mat)
            mouse(imouse).task(itask).pos(ipos).ind = [];
            for iexp = 1:length(mouse(imouse).task(itask).ind)
                if xd.Position(mouse(imouse).task(itask).ind(iexp)) == pv.pos_mat(ipos);
                    mouse(imouse).task(itask).pos(ipos).ind = [mouse(imouse).task(itask).pos(ipos).ind mouse(imouse).task(itask).ind(iexp)];
                end
            end
        end
    end
    
    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).task(itask).block(1).half(1).nCorrect = 0;
        mouse(imouse).task(itask).block(1).half(2).nCorrect = 0;
        mouse(imouse).task(itask).block(2).half(1).nCorrect = 0;
        mouse(imouse).task(itask).block(2).half(2).nCorrect = 0;
        mouse(imouse).task(itask).block(1).half(1).nMissed = 0;
        mouse(imouse).task(itask).block(1).half(2).nMissed = 0;
        mouse(imouse).task(itask).block(2).half(1).nMissed = 0;
        mouse(imouse).task(itask).block(2).half(2).nMissed = 0;
        for ipos = 1:length(pv.pos_mat)
            iexp = 1;
            mouse(imouse).task(itask).pos(ipos).block(1).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(1).nMissed = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).nMissed = 0;
            mouse(imouse).task(itask).pos(ipos).block(1).half(1).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(1).half(1).nMissed = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).half(1).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).half(1).nMissed = 0;
            mouse(imouse).task(itask).pos(ipos).block(1).half(2).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(1).half(2).nMissed = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).half(2).nCorrect = 0;
            mouse(imouse).task(itask).pos(ipos).block(2).half(2).nMissed = 0;
            for ind = mouse(imouse).task(itask).pos(ipos).ind
                fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(ind)) '-' cell2mat(xd.DateStr(ind)) '.mat']);
                load(fn)
                trials = cell2mat(xd.TrialRangeToUse(ind));
                if isempty(trials)
                    trials = input.trialSinceReset;
                    range = 1:trials;
                else
                    range = str2num(trials);
                end
                midpoint = ceil((range(end)-range(1))/2);
                firsthalftrials = range(1):midpoint;
                secondhalftrials = midpoint+1:range(end);
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==0), range);
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==1), range);
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).intensities = unique(chop(unique(cell2mat(input.gratingContrast)),4));                
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities = unique(chop(unique(cell2mat(input.gratingContrast(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).trial_ind))),4));                    
                fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(ind)) '-max.mat']);
                load(fn);
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).thresh = fitS(1).thresh;
                mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).ci95 = fitS(1).bootStats.ci95;
                ci95_ratio_1 = fitS(1).bootStats.ci95(2)./fitS(1).bootStats.ci95(1);
                if ci95_ratio_1 < 2
                    mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).well_fit = 1;
                else
                    mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).well_fit = 0;
                end
                for iblock = 1:2
                    for icon = 1:length(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities)
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nTrials(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).trial_ind, find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon))));
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nCorrect(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon)), intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nMissed(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon)), find(cell2mat(input.reactTimesMs)>550))));
                        [x y] = binofit(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nCorrect(icon),mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nCorrect(icon)+mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nMissed(icon));
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).pctCorrect = x;
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).ci95 = y;
                        for ihalf = 1:2
                            if ihalf  == 1
                                half_trials = intersect(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).trial_ind,firsthalftrials);
                            elseif ihalf  == 2
                                half_trials = intersect(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).trial_ind,secondhalftrials);
                            end
                            mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nTrials(icon) = length(intersect(half_trials, find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon))));
                            mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nCorrect(icon) = length(intersect(half_trials, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon)), intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                            mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nMissed(icon) = length(intersect(half_trials, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities(icon)), find(cell2mat(input.reactTimesMs)>550))));
                            [x y] = binofit(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nCorrect(icon),mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nCorrect(icon)+mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nMissed(icon));
                            mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).half(ihalf).pctCorrect = x;
                            mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).binofit(icon).half(ihalf).ci95 = y;
                        end
                    end
                    if ci95_ratio_1 < 2
                        mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect = mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nCorrect(1);
                        mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed = mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).nMissed(1);
                        for ihalf = 1:2
                            mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect = mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nCorrect(1);
                            mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed = mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(iblock).half(ihalf).nMissed(1);
                        end
                    end
                end
                for icon = 1:length(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).intensities)
                    mouse(imouse).task(itask).pos(ipos).expt(iexp).p(icon) = binocdf(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).nCorrect(icon), (mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).nCorrect(icon) + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).nMissed(icon)),  mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).binofit(icon).pctCorrect);                            
                    for ihalf = 1:2
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).half(ihalf).p(icon) = binocdf(mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).half(ihalf).nCorrect(icon), (mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).half(ihalf).nCorrect(icon) + mouse(imouse).task(itask).pos(ipos).expt(iexp).block(2).half(ihalf).nMissed(icon)),  mouse(imouse).task(itask).pos(ipos).expt(iexp).block(1).binofit(icon).half(ihalf).pctCorrect);                            
                    end
                end
                iexp = 1+iexp;
            end
            for iblock = 1:2
                [x y] = binofit(mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect,mouse(imouse).task(itask).pos(ipos).block(iblock).nCorrect+mouse(imouse).task(itask).pos(ipos).block(iblock).nMissed);
                mouse(imouse).task(itask).pos(ipos).block(iblock).pctCorrect = x;
                mouse(imouse).task(itask).pos(ipos).block(iblock).ci95 = y;
                for ihalf = 1:2
                    [x y] = binofit(mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect,mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nCorrect+mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).nMissed);
                    mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).pctCorrect = x;
                    mouse(imouse).task(itask).pos(ipos).block(iblock).half(ihalf).ci95 = y;
                end
            end
        end
        mouse(imouse).task(itask).pos(1).p = binocdf(mouse(imouse).task(itask).pos(1).block(2).nCorrect, (mouse(imouse).task(itask).pos(1).block(2).nCorrect + mouse(imouse).task(itask).pos(1).block(2).nMissed),  mouse(imouse).task(itask).pos(2).block(1).pctCorrect);                            
        mouse(imouse).task(itask).pos(2).p = binocdf(mouse(imouse).task(itask).pos(2).block(1).nCorrect, (mouse(imouse).task(itask).pos(2).block(1).nCorrect + mouse(imouse).task(itask).pos(2).block(1).nMissed),  mouse(imouse).task(itask).pos(1).block(2).pctCorrect);                                
        for ihalf = 1:2
            mouse(imouse).task(itask).pos(1).half(ihalf).p = binocdf(mouse(imouse).task(itask).pos(1).block(2).half(ihalf).nCorrect, (mouse(imouse).task(itask).pos(1).block(2).half(ihalf).nCorrect + mouse(imouse).task(itask).pos(1).block(2).half(ihalf).nMissed),  mouse(imouse).task(itask).pos(2).block(1).half(ihalf).pctCorrect);                            
            mouse(imouse).task(itask).pos(2).half(ihalf).p = binocdf(mouse(imouse).task(itask).pos(2).block(1).half(ihalf).nCorrect, (mouse(imouse).task(itask).pos(2).block(1).half(ihalf).nCorrect + mouse(imouse).task(itask).pos(2).block(1).half(ihalf).nMissed),  mouse(imouse).task(itask).pos(1).block(2).half(ihalf).pctCorrect);                                
        end
        for iblock = 1:2
            for ihalf = 1:2
                mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect = mouse(imouse).task(itask).pos(1).block(iblock).half(ihalf).nCorrect +mouse(imouse).task(itask).pos(2).block(iblock).half(ihalf).nCorrect;
                mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed = mouse(imouse).task(itask).pos(1).block(iblock).half(ihalf).nMissed +mouse(imouse).task(itask).pos(2).block(iblock).half(ihalf).nMissed;
                mouse(imouse).task(itask).block(iblock).half(ihalf).pctCorrect = mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect/(mouse(imouse).task(itask).block(iblock).half(ihalf).nMissed+mouse(imouse).task(itask).block(iblock).half(ihalf).nCorrect);
            end
        end
    end
end
    
    
                
    