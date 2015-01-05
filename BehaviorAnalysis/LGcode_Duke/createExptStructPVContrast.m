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
        mouse(imouse).name = pv.mouse_mat(imouse);
        mouse(imouse).task(itask).ind = intersect(find(strcmp(xd.Task,task)), mouse(imouse).ind);
        for ipos = 1:length(pv.pos_mat)
            pos = pv.pos_mat(ipos);
            mouse(imouse).task(itask).pos(ipos).ind = intersect(find(xd.Position==pos), mouse(imouse).task(itask).ind);
            mouse(imouse).task(itask).pos(ipos).name = pos;
        end
    end
    
    % find all days for each power
    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                pow = pv.power_mat(ipow);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).ind = intersect(find(xd.Power==pow), mouse(imouse).task(itask).pos(ipos).ind);
                mouse(imouse).task(itask).pos(ipos).pow(ipow).name = pow;
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
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).name = con;
                end
            end
        end    
    end


    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                for icon = 1:length(pv.basecon_mat)
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt = [];
                    nexp = length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind);
                    thresh = zeros(nexp,2);
                    early_pct = zeros(nexp,2);
                    max = zeros(nexp,2);
                    fit_early_max = zeros(nexp,3,2);
                    date = zeros(nexp,1);
                    for iexp = 1:nexp
                        dateStr = cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp)));
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).date = dateStr;
                        fn = fullfile(rc.fitOutputMatDir, ['subj' num2str(pv.mouse_mat(imouse)) '-' dateStr '-' cell2mat(xd.DataBlock(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '.mat']);
                        load(fn);
                        for las = 1:length(fitS)
                            %check confidence intervals
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).ci95_ratio = fitS(las).bootStats.ci95(2)./fitS(las).bootStats.ci95(1);
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).ci95_ratio < 2
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 1;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,1) = 1;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,1) = 0;
                            end
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh = fitS(las).thresh;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max = fitS(las).coefEsts(3);
                            %check lapse rate
                            if fitS(las).coefEsts(3)<.9
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,2) = 0;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,2) = 1;
                            end
                            %check false alarm rate
                            ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp);
                            fName = rc.computeFName(pv.mouse_mat(imouse), cell2mat(xd.DateStr(ind)));
                            n  = dir([fName(1:(length(rc.pathStr)+17)) '*']);
                            clear ds;
                            clear temp;
                            if size(n,1) == 1
                                fName = fullfile(rc.pathStr, n.name);
                                ds = mwLoadData(fName);
                            elseif size(n,1) > 1
                                if xd.MergeAllMatFiles(ind) == 1
                                    for ifile = 1:size(n,1)
                                        fName = fullfile(rc.pathStr, n(ifile).name);
                                        temp = mwLoadData(fName);
                                        ds(ifile) = temp;
                                        if ifile>1
                                            ds = concatenateDataBlocks(ds);
                                        end
                                    end
%                                 elseif ~isnan(xd.ChooseMatFile(ind))
%                                     for ifile = 1:size(n,1)
%                                         if ~isempty(strfind(n(ifile).name,num2str(xd.ChooseMatFile(ind))))
%                                             fName = fullfile(rc.pathStr, n(ifile).name);
%                                             ds = load(fName);
%                                         end
%                                     end
                                else
                                    error('Too many mat files- need to choose');
                                end
                            end    
                            trials = str2num(cell2mat(xd.TrialRangeToUse(ind)));
                            if ~isnan(trials)
                                ds = trialChopper(ds,[trials(1) trials(end)]);
                            end
                            ind_all = find(cell2mat(ds.tBlock2TrialNumber)==(las-1));
                            ind_hold = intersect(find(cell2mat(ds.holdTimesMs)>200), ind_all);
                            RTs = cell2mat(ds.reactTimesMs);
                            ind_early = find(RTs(ind_hold)<=125);
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct = length(ind_early)./length(ind_hold);
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct>0.5
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good = 0;
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,3) = 0;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,3) = 1;
                            end
                            fit_early_max(iexp,:,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max;
                            %concatenate all good days                           
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good
                               thresh(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh;
                               early_pct(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct;
                               max(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max;
                            else
                               thresh(iexp,las) = NaN;
                               early_pct(iexp,las) = NaN;
                               max(iexp,las) = NaN;
                            end
                        end
                    end
                    good_ind = find(sum(~isnan(thresh),2)==2);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_ind = good_ind;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh = thresh;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).early_pct = early_pct;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).max = max;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max = fit_early_max;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).date = xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind);
                    if length(good_ind)>0
                        figure;
                        subplot(2,3,1)
                        plot([1:2],thresh,'o-')
                        ylim([0 0.5])
                        title('thresh')
                        subplot(2,3,4)
                        plot(1:nexp,thresh,'o')
                        ylim([0 0.5])
                        xlim([1 nexp])
                        title('Threshold- control(blue) LED(green)')
                        subplot(2,3,2)
                        plot([1:2],early_pct,'-')
                        ylim([0 0.5])
                        title('early pct')
                        subplot(2,3,3)
                        plot([1:2],max,'-')
                        ylim([0.85 1.0])
                        title('max')
                        subplot(2,3,5)
                        plot(1:nexp,(thresh(:,2)./thresh(:,1))','o')
                        ylim([0.5 3])
                        xlim([1 nexp])
                        title('Threshold- LED/control')
                        subplot(2,3,6)
                        plot((thresh(:,2)./thresh(:,1)),(early_pct(:,2)./early_pct(:,1)),'o')
                        ylim([0.5 1.5])
                        xlim([0.5 3])
                        title('LED/Control: Threshold(X) vs FalseAlarm (Y)')
                        suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
                        fname = ['ThresholdSummary_i' num2str(pv.mouse_mat(imouse)) '_' pv.task_mat(itask,:) num2str(pv.power_mat(ipow)) 'mW.pdf'];
                        fn_out = fullfile(rc.fitOutputSummary, fname);
                        exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                    end
                    if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
                        figure;
                        subplot(2,3,4); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1,2),'o'); title('LED')
                        subplot(2,3,1); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1,1),'o'); title('fit')
                        subplot(2,3,2); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2,1),'o'); title('early')
                        subplot(2,3,5); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2,2),'o'); title('LED')
                        subplot(2,3,3); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,3,1),'o'); title('max')
                        subplot(2,3,6); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,3,2),'o'); title('LED')
                        suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
                        fn_out = fullfile(rc.fitOutputSummary, ['FitEarlyMax_i' num2str(pv.mouse_mat(imouse)) '_'  pv.task_mat(itask,:) '_' num2str(pv.power_mat(ipow)) 'mW.pdf']);
                        exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                        figure;
                        for las = 1:2
                            fit_early = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1:2,las),2))==2);
                            fit_max = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,[1 3],las),2))==2);
                            early_max = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2:3,las),2))==2);
                            sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1:2,las),2)
                            size(fit_early)
                            size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1)
                            subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1),fit_early,'o'); title('fit_early')
                            subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), fit_max,'o'); title('fit_max')
                            subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), early_max,'o'); title('early_max')
                        end
                        suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
                        fn_out = fullfile(rc.fitOutputSummary, ['FitEarlyMax_pairs_i' num2str(pv.mouse_mat(imouse)) '_'  pv.task_mat(itask,:) '_' num2str(pv.power_mat(ipow)) 'mW.pdf']);
                        exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                    end
                end
            end
        end
    end
end
fn_out = fullfile(rc.fitOutputSummary, ['PVcon_summarystruct_' date '.mat']);
save(fn_out, 'mouse');