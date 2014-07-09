function mouse = createExptStructPV(tasks)
%takes all data from xls, sorts according to task type and gets thresholds
%etc

rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

pv = behavParamsPV;
for itask = tasks
    for imouse = 1:length(pv.mouse_mat)
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

    for imouse = 1:length(pv.mouse_mat)
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

    if itask > 1;
        for imouse = 1:length(pv.mouse_mat)
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

    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt = [];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(1).nMissed = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).nCorrect = [0 0];
                mouse(imouse).task(itask).pos(ipos).pow(ipow).block(2).nMissed = [0 0];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
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
                        if fitS(1).coefEsts(3)>=.9
                            if fitS(2).coefEsts(3)>=.9
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
                            end
                        else
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
                        end
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1;
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                                    fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                                    load(fn)
                                    trials = cell2mat(xd.TrialRangeToUse(ind));
                                    if isempty(trials)
                                        trials = input.trialSinceReset;
                                        range = 1:trials;
                                    else
                                        range = str2num(trials);
                                    end
                                    block1_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)==0), range);
                                    block2_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)>0), range);
                                    block1_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block1_ind_all);
                                    block2_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block2_ind_all);
                                    RTs = cell2mat(input.reactTimesMs);
                                    block1_ind_early = find(RTs(block1_ind_good)<=125);
                                    block2_ind_early = find(RTs(block2_ind_good)<=125);
                                    early_pct = sum([length(block1_ind_early) length(block2_ind_early)])./length(range);
                                    if early_pct<0.5
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 1;
                                    else
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 0;
                                    end
                                end
                            end
                        end
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
                            if fitS(1).coefEsts(3)>=.9
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
                            end
                        end
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                                fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                                load(fn)
                                trials = cell2mat(xd.TrialRangeToUse(ind));
                                if isempty(trials)
                                    trials = input.trialSinceReset;
                                    range = 1:trials;
                                else
                                    range = str2num(trials);
                                end
                                block1_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)==0), range);
                                block2_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)>0), range);
                                block1_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block1_ind_all);
                                block2_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block2_ind_all);
                                RTs = cell2mat(input.reactTimesMs);
                                block1_ind_early = find(RTs(block1_ind_good)<=125);
                                block2_ind_early = find(RTs(block2_ind_good)<=125);
                                early_pct = sum([length(block1_ind_early) length(block2_ind_early)])./length(range);
                                if early_pct<0.5
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 1;
                                else
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 0;
                                end
                            end
                        end
                    end
                end     
%                 if itask > 1
%                     for icon = 1:length(pv.basecon_mat)
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt = [];
%                         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                             mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser = [];
%                             fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '.mat']);
%                             load(fn);
%                             if length(fitS) == 2
%                                 for las = 1:length(fitS)
%                                     ci95_ratio(las) = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
%                                     if ci95_ratio(las) < 2
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 1;
%                                     else
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
%                                     end
%                                     mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh = fitS(las).thresh;
%                                     mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
%                                 end
%                                 if fitS(1).coefEsts(3)>=.9
%                                     if fitS(2).coefEsts(3)>=.9
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 1;
%                                     else
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
%                                     end
%                                 else
%                                     mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
%                                 end
%                                 if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                                     if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1;
%                                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
%                                             fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
%                                             load(fn)
%                                             trials = cell2mat(xd.TrialRangeToUse(ind));
%                                             if isempty(trials)
%                                                 trials = input.trialSinceReset;
%                                                 range = 1:trials;
%                                             else
%                                                 range = str2num(trials);
%                                             end
%                                             block1_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)==0), range);
%                                             block2_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)>0), range);
%                                             block1_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block1_ind_all);
%                                             block2_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block2_ind_all);
%                                             RTs = cell2mat(input.reactTimesMs);
%                                             block1_ind_early = find(RTs(block1_ind_good)<=125);
%                                             block2_ind_early = find(RTs(block2_ind_good)<=125);
%                                             early_pct = sum([length(block1_ind_early) length(block2_ind_early)])./length(range);
%                                             if early_pct<0.5
%                                                 mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 1;
%                                             else
%                                                 mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 0;
%                                             end
%                                         end
%                                     end
%                                 end
%                             elseif length(fitS) == 1
%                                 mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh = NaN;
%                                 mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).max = NaN;
%                                 mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good = 0;
%                                 for las = 1:length(fitS)
%                                     ci95_ratio(las) = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
%                                     mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh = fitS(las).thresh;
%                                     mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
%                                     if ci95_ratio(las) < 2
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 1;
%                                     else
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
%                                     end
%                                     if fitS(1).coefEsts(3)>=.9
%                                     	mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 1;
%                                     else
%                                         mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss = 0;
%                                     end
%                                 end
%                                 if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                                     if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
%                                         fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
%                                         load(fn)
%                                         trials = cell2mat(xd.TrialRangeToUse(ind));
%                                         if isempty(trials)
%                                             trials = input.trialSinceReset;
%                                             range = 1:trials;
%                                         else
%                                             range = str2num(trials);
%                                         end
%                                         block1_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)==0), range);
%                                         block2_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)>0), range);
%                                         block1_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block1_ind_all);
%                                         block2_ind_good = intersect(find(cell2mat(input.holdTimesMs)>200), block2_ind_all);
%                                         RTs = cell2mat(input.reactTimesMs);
%                                         block1_ind_early = find(RTs(block1_ind_good)<=125);
%                                         block2_ind_early = find(RTs(block2_ind_good)<=125);
%                                         early_pct = sum([length(block1_ind_early) length(block2_ind_early)])./length(range);
%                                         if early_pct<0.5
%                                             mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 1;
%                                         else
%                                             mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = 0;
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
            end
        end
    end
end