function mouse = createExptStructAll
%takes all data from xls, sorts according to task type and gets thresholds
%etc

    rc = behavConstsHADC8;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAll;

    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        for itask = 1:length(pv.task_mat)
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
    end

    for imouse = 1:length(pv.mouse_mat)
        for itask = 1:length(pv.task_mat)
            for ipos = 1:length(pv.pos_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
                    fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).ind(iexp))) '.mat']);
                    load(fn);    
                    ci95_ratio_1 = fitS(1).bootStats.ci95(2)./fitS(1).bootStats.ci95(1);
                    if ci95_ratio_1 < 2
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh = fitS(1).thresh;
                    else
                        mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh = NaN;
                    end
                end
            end
        end
    end
    
    for imouse = 1:length(pv.mouse_mat)
        for itask = pv.task_mat
            for ipos = 1:length(pv.pos_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
                    fn = fullfile('/Users/lindsey/Desktop/Data', ['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).ind(iexp))) '.mat']);
                    load(fn);
                    trials = cell2mat(xd.TrialRangeToUse(mouse(imouse).task(itask).pos(ipos).ind(iexp)));
                    if isempty(trials)
                        trials = input.trialSinceReset;
                        range = 1:trials;
                    else
                        range = str2num(trials);
                    end
                    block1_ind_all = intersect(find(cell2mat(input.tTrialLaserPowerMw)==0),range);
                    RTs = cell2mat(input.reactTimesMs);
                    block1_ind_early = find(RTs(block1_ind_all)<=125);
                    mouse(imouse).task(itask).pos(ipos).expt(iexp).n_all = length(block1_ind_all);
                    mouse(imouse).task(itask).pos(ipos).expt(iexp).n_early = length(block1_ind_early);   
                end
            end
        end
    end
end