function mouse = plotProbeDataByHalf(tasks)

rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

pv = behavParamsPV;

for itask = tasks
    for imouse = [2 4]
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

    for imouse = [2 4]
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



    for imouse = [2 4]
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt = [];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).half(1).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).half(1).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).half(1).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).half(1).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).half(2).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).half(2).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).half(2).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).half(2).nMissed = [0 0];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser = [];
                    fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                    load(fn);
                    if length(fitS) == 1
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).b2probe = 1;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh = NaN;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).max = NaN;
                        ci95_ratio_1 = fitS(1).bootStats.ci95(2)./fitS(1).bootStats.ci95(1);
                        if ci95_ratio_1 < 2
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh = fitS(1).thresh;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).max = fitS(1).coefEsts(3);
                            fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                            load(fn)
                            trials = cell2mat(xd.TrialRangeToUse(ind));
                            if isempty(trials)
                                trials = input.trialSinceReset;
                                range = 1:trials;
                            else
                                range = str2num(trials);
                            end
                            firsthalftrials = find(cell2mat(input.tTotalReqHoldTimeMs)<=1000);
                            secondhalftrials = find(cell2mat(input.tTotalReqHoldTimeMs)>1000);
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(1).half(1).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==0), intersect(firsthalftrials,range));
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(1).half(2).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==0), intersect(secondhalftrials,range));
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).half(1).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==1), intersect(firsthalftrials,range));
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).half(2).trial_ind = intersect(find(cell2mat(input.tBlock2TrialNumber) ==1), intersect(secondhalftrials,range));
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(1).intensities = unique(chop(unique(cell2mat(input.gratingContrast)),4));                
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities = unique(chop(unique(cell2mat(input.gratingContrast(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).half(1).trial_ind))),4));                    
                            for iblock = 1:2
                                for ihalf = 1:2
                                    for icon = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities)
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nTrials(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).trial_ind, find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon))));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nCorrect(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon)), intersect(find(cell2mat(input.reactTimesMs)>125), find(cell2mat(input.reactTimesMs)<550)))));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nMissed(icon) = length(intersect(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).trial_ind, intersect(find(chop(cell2mat(input.gratingContrast),4)==mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(2).intensities(icon)), find(cell2mat(input.reactTimesMs)>550))));
                                        [x y] = binofit(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nCorrect(icon),mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nCorrect(icon)+mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nMissed(icon));
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).binofit(icon).pctCorrect = x;
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).binofit(icon).ci95 = y;
                                    end
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nCorrect = mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nCorrect + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nCorrect;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nMissed = mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nMissed + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).block(iblock).half(ihalf).nMissed;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).b2probe = 1;
                                end
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

for imouse = [2 4]
    for ipos = 1:length(pv.pos_mat)
        if mouse(imouse).task(itask).pos(ipos).ind>0
            figure;
            start = 1;
            for ipow = 1:length(pv.power_mat)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).ind>0
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).b2probe == 1
                        subplot(2,2,start)
                        col_str = strvcat('k','c');
                        for iblock = 1:2
                            for ihalf = 1:2
                                pctCorrect = mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nCorrect(1)./(mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nCorrect(1)+mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nMissed(1));
                                plot(ihalf, pctCorrect, ['s' col_str(iblock,:)])
                                hold on
                                text(ihalf, iblock*.1, num2str(mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nCorrect(1)+mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).half(ihalf).nMissed(1)), 'Color', col_str(iblock))
                                hold on
                            end
                        end
                        xlim([0 3])
                        ylim([0 1])
                        title(['Power = ' num2str(pv.power_mat(ipow))])
                        start = start+1;
                    end
                end
            end
            suptitle(['Mouse i' num2str(pv.mouse_mat(imouse)) '  Pos = ' num2str(pv.pos_mat(ipos))])
        end
    end
end
                        