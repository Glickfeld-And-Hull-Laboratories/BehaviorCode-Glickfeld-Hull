function mouse = createExptStructArch(tasks)
%takes all data from xls, sorts according to task type and gets thresholds
%etc

rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

pv = behavParamsPV;
for itask = tasks
    for imouse = find(pv.arch_mat == 1)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).task(itask).ind = intersect(find(xd.Task == itask), mouse(imouse).ind);
        for ipos = 1:length(pv.pos_mat)
            mouse(imouse).task(itask).pos(ipos).ind = [];
            for iexp = 1:length(mouse(imouse).task(itask).ind)
                if str2num(cell2mat(xd.Position(mouse(imouse).task(itask).ind(iexp)))) == pv.pos_mat(ipos);
                    mouse(imouse).task(itask).pos(ipos).ind = [mouse(imouse).task(itask).pos(ipos).ind mouse(imouse).task(itask).ind(iexp)];
                end
            end
        end
    end

    for imouse = find(pv.arch_mat == 1)
        for ipos = 1:length(pv.pos_mat)
            ind = mouse(imouse).task(itask).pos(ipos).ind;
            for ipow = 1:length(pv.power_mat)
                mouse(imouse).task(itask).pos(ipos).pow(ipow).ind = [];
                for iexp = 1:length(ind)
                    if xd.Power(ind(iexp)) == pv.power_mat(ipow);
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).ind = [mouse(imouse).task(itask).pos(ipos).pow(ipow).ind ind(iexp)];
                        end
                    end
                end
            end
        end    
    end

    if itask > 1;
        for imouse = find(pv.arch_mat == 1)
            for ipos = 1:length(pv.pos_mat)
                for ipow = 1:length(pv.power_mat)
                    ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind;
                    for icon = 1:length(pv.basecon_mat);
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind = [];
                        for iexp = 1:length(ind)
                            if xd.BaseContrast(ind(iexp)) == pv.basecon_mat(icon);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind = [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind ind(iexp)];
                            end
                        end
                    end
                end
            end    
        end
    end

    for imouse = find(pv.arch_mat == 1)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).b2probe = 0;
                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt = [];
                ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind;
                for iexp = 1:length(ind)
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).AC = xd.SplitBlock1(ind(iexp));
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser = [];
                    fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                    load(fn);
                    if length(fitS) == 2
                        for las = 1:length(fitS)
                            ci95_ratio(las) = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
                            if ci95_ratio(las) < 2
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 0;
                            end
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).thresh = fitS(las).thresh;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
                        end
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).b2probe = 0;
                    elseif length(fitS) == 1
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh = NaN;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).max = NaN;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good = 0;
                        for las = 1:length(fitS)
                            ci95_ratio(las) = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).thresh = fitS(las).thresh;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
                            if ci95_ratio(las) < 2
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 0;
                            end
                            fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                            load(fn)
                            if itask == 1
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).b2probe = 1;
                                trials = cell2mat(xd.TrialRangeToUse(ind(iexp)));
                                if isempty(trials)
                                    trials = input.trialSinceReset;
                                    range = 1:trials;
                                else
                                    range = str2num(trials);
                                end
                                midpoint = ceil((range(end)-range(1))/2);
                                firsthalftrials = range(1):midpoint;
                                secondhalftrials = midpoint+1:range(end);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(1).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==0), range);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==1), range);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(1).intensities = unique(chop(unique(cell2mat(input.gratingContrast)),4));                
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities = unique(chop(unique(cell2mat(input.gratingContrast(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).trial_ind))),4));                    
                                for iblock = 1:2
                                    for icon = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities)
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nTrials(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).trial_ind, find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon))));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nCorrect(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon)), intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nMissed(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon)), find(cell2mat(input.reactTimesMs)>550))));
                                        [x y] = binofit(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nCorrect(icon),mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nCorrect(icon)+mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nMissed(icon));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).binofit(icon).pctCorrect = x;
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).binofit(icon).ci95 = y;
                                    end
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nCorrect = mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nCorrect + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nCorrect;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nMissed = mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nMissed + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).nMissed;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).b2probe = 1;
                                end
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh = NaN;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).max = NaN;
                            end
                        end
                    end                
                end
            end
        end
    end
end