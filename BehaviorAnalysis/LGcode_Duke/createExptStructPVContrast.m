function mouse = createExptStructPVContrast(tasks)
%takes all data from xls, sorts according to task type and gets thresholds
%etc

rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

pv = behavParamsPVContrast;
for itask = tasks
    task = pv.task_mat(itask,:);
    % find all days for each position
    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).task(itask).ind = intersect(find(xd.Task == task), mouse(imouse).ind);
        for ipos = 1:length(pv.pos_mat)
            mouse(imouse).task(itask).pos(ipos).ind = [];
            for iexp = 1:length(mouse(imouse).task(itask).ind)
                if str2num(cell2mat(xd.Position(mouse(imouse).task(itask).ind(iexp)))) == pv.pos_mat(ipos);
                    mouse(imouse).task(itask).pos(ipos).ind = [mouse(imouse).task(itask).pos(ipos).ind mouse(imouse).task(itask).ind(iexp)];
                end
            end
        end
    end
    
    % find all days for each power
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

    % find all days for each base contrast
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


    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt = [];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser = [];
                    fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                    load(fn);
                    for las = 1:length(fitS)
                        %check confidence intervals
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).ci95_ratio = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).ci95_ratio < 2
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 1;
                        else
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 0;
                        end
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).thresh = fitS(las).thresh;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
                        %check lapse rate
                        if fitS(las).coefEsts(3)>=.9
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 1;
                        else
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 0;
                        end
                        %check false alarm rate
                        fn = fullfile(rc.pathStr, ['data-i' num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                        load(fn)
                        trials = cell2mat(xd.TrialRangeToUse(ind));
                        if ~isempty(trials)
                            input = trialChopper(input,trials);
                        end
                        ind_all = find(cell2mat(input.tBlock2TrNumber==(las-1)));
                        ind_hold = intersect(find(cell2mat(input.holdTimesMs)>200), ind_all);
                        RTs = cell2mat(input.reactTimesMs);
                        ind_early = find(RTs(ind_hold)<=125);
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).early_pct = length(ind_early)./length(ind_hold);
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).early_pct<0.5
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 1;
                        else
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(las).good = 0;
                        end
                    end
                end
            end
        end
    end
end