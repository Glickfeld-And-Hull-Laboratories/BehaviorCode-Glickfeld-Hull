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
        mouse(imouse).task(itask).ind = intersect(find(strcmp(xd.Task,task)), mouse(imouse).ind);
        for ipos = 1:length(pv.pos_mat)
            pos = pv.pos_mat(ipos);
            mouse(imouse).task(itask).pos(ipos).ind = intersect(find(xd.Position==pos), mouse(imouse).task(itask).ind);
        end
    end
    
    % find all days for each power
    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                pow = pv.power_mat(ipow);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).ind = intersect(find(xd.Power==pow), mouse(imouse).task(itask).pos(ipos).ind);
            end
        end    
    end

    % find all days for each base contrast
    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                for icon = 1:length(pv.basecon_mat);
                    con = pv.basecon_mat(icon);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind = intersect(find(xd.BaseContrast==con), mouse(imouse).task(itask).pos(ipos).pow(ipow).ind);
                end
            end
        end    
    end


    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                for icon = 1:length(pv.basecon_mat)
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt = [];
                    thresh = [];
                    early_pct = [];
                    max = [];
                    days = [];
                    for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser = [];
                        fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '.mat']);
                        load(fn);
                        for las = 1:length(fitS)
                            %check confidence intervals
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).ci95_ratio = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).ci95_ratio < 2
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                            end
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh = fitS(las).thresh;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
                            %check lapse rate
                            if fitS(las).coefEsts(3)<.9
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                            end
                            %check false alarm rate
                            ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp);
                            fName = rc.computeFName(pv.mouse_mat(imouse), cell2mat(xd.DateStr(ind)));
                            n  = dir([fName(1:(length(rc.pathStr)+17)) '*']);
                            if size(n,1) == 1
                                fName = fullfile(rc.pathStr, n.name);
                                ds = load(fName);
                            elseif size(n,1) > 1
                                if uo.MergeMats == 1
                                    for ifile = 1:size(n,1)
                                        fName = fullfile(rtc.pathStr, n(ifile).name);
                                        ds = load(fName);
                                        ds = concatenateDataBlocks(ds);
                                    end
                                else
                                    error('Too many mat files- need to choose');
                                end              
                            end    
                            trials = str2num(cell2mat(xd.TrialRangeToUse(ind)));
                            if ~isnan(trials)
                                ds.input = trialChopper(ds.input,[trials(1) trials(end)]);
                            end
                            ind_all = find(cell2mat(ds.input.tBlock2TrialNumber)==(las-1));
                            ind_hold = intersect(find(cell2mat(ds.input.holdTimesMs)>200), ind_all);
                            RTs = cell2mat(ds.input.reactTimesMs);
                            ind_early = find(RTs(ind_hold)<=125);
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct = length(ind_early)./length(ind_hold);
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct>0.5
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                            end
                        end
                        %concatenate all good days
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good
                                thresh = [thresh; mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh];
                                early_pct = [early_pct; mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).early_pct  mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).early_pct];
                                max = [max; mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).max mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).max];
                                days = [days; iexp];
                            end
                        end
                    end
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).n = length(thresh);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_ind = days;
                    if length(thresh)>0
                        figure;
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_all = thresh;
                        subplot(2,3,1)
                        plot([1:2],thresh,'-')
                        title('thresh')
                        subplot(2,3,4)
                        plot(days,thresh','o')
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).n = length(thresh);
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).early_pct_all = early_pct;
                        subplot(2,3,2)
                        plot([1:2],early_pct,'-')
                        title('early pct')
                        subplot(2,3,5)
                        plot(days,early_pct','o')
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).max_all = max;
                        subplot(2,3,3)
                        plot([1:2],max,'-')
                        title('max')
                        subplot(2,3,6)
                        plot(days,max','o')
                        suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' pv.power_mat(ipow)]);
                    end
                end
            end
        end
    end
end

%% combining multiple mat files

function outS = concatenateDataBlocks(blockC)
    fNames = fieldnames(blockC(1));
    nF = length(fNames);
    outS = struct([]);
    trs = blockC(1).trialSinceReset;
    for iF = 1:nF
        fN = cell2mat(fNames(iF,:));
        outS(1).(fN) = blockC(1).(fN);
        if size(blockC(1).(fN),2) == trs;
            for iblock = 2:size(blockC,2)
                outS.(fN) = [outS.(fN) blockC(iblock).(fN)];
            end
        else
            outS.(fN) = blockC(1).(fN);
        end
    end