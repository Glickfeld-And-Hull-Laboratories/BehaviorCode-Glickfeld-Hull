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
                    good_thresh = zeros(nexp,2);
                    good_early_pct = zeros(nexp,2);
                    good_max = zeros(nexp,2);
                    fit_early_max = zeros(nexp,3,2);
                    date = zeros(nexp,1);
                    for iexp = 1:nexp
                        dateStr = cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp)));
                        mouseStr = num2str(xd.Subject(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp)));
                        disp([dateStr ' i' mouseStr])
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
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,3) = 0;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,3) = 1;
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
                                elseif xd.ChooseMatFile(ind)>1
                                    for ifile = 1:size(n,1)
                                        if ~isempty(strfind(n(ifile).name,num2str(xd.ChooseMatFile(ind))))
                                            fName = fullfile(rc.pathStr, n(ifile).name);
                                            ds = mwLoadData(fName);
                                        end
                                    end
                                else
                                    doLever = zeros(size(n,1));
                                    for ifile = 1:size(n,1)
                                        fName = fullfile(rc.pathStr, n(ifile).name);
                                        input = mwLoadData(fName);
                                        doLever(ifile,1) = input.doLever;
                                    end
                                    doLever_list = find(doLever);
                                    if size(doLever_list,1) == 1
                                        fName = fullfile(rc.pathStr, n(doLever_list).name);
                                        ds = mwLoadData(fName);
                                    else
                                        error('Too many mat files- need to choose');
                                    end
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
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,2) = 0;
                            else
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max(1,2) = 1;
                            end
                            fit_early_max(iexp,:,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).fit_early_max;
                            %concatenate all days
                            thresh(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh;
                            early_pct(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct;
                            max(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max;
                            %concatenate all good days                           
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).good
                               good_thresh(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).thresh;
                               good_early_pct(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).early_pct;
                               good_max(iexp,las) = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(las).max;
                            else
                               good_thresh(iexp,las) = NaN;
                               good_early_pct(iexp,las) = NaN;
                               good_max(iexp,las) = NaN;
                            end
                        end
                    end
                    both_good_ind = find(sum(~isnan(good_thresh),2)==2);
                    good_ind = find(sum(~isnan(good_thresh),2)>=1);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_ind = good_ind;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).both_good_ind = both_good_ind;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_thresh = good_thresh;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_early_pct = good_early_pct;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).good_max = good_max;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh = thresh;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).early_pct = early_pct;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).max = max;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max = fit_early_max;
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).date = xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind);
                    if length(good_ind)>0
                        figure;
                        subplot(3,3,1)
                        plot([1:2],good_thresh,'o-')
                        title('thresh')
                        subplot(3,3,4)
                        plot(good_ind,good_thresh(good_ind,:),'o')
                        title('Control(blue) LED(green)')
                        ylabel('Threshold')
                        subplot(3,3,2)
                        plot([1:2],early_pct,'-')
                        ylim([0 1.0])
                        title('early pct')
                        subplot(3,3,3)
                        plot([1:2],max,'-')
                        ylim([0 1.0])
                        title('max')
                        subplot(3,3,5)
                        plot(good_ind,(good_thresh(good_ind,2)./good_thresh(good_ind,1))','o')
                        ylim([0.5 3])
                        ylabel('Threshold Ratio')
                        title('LED/Control')
                        subplot(3,3,6)
                        plot((good_early_pct(:,2)./good_early_pct(:,1)),(good_thresh(:,2)./good_thresh(:,1)),'o')
                        xlim([0.5 1.5])
                        ylim([0.5 3])
                        xlabel('FArate Ratio')
                        ylabel('Threshold Ratio')
                        title('LED/Control')
                        subplot(3,3,7)
                        plot(mean(early_pct,2),(good_thresh(:,2)./good_thresh(:,1)),'o')
                        xlabel('FArate')
                        ylabel('Threshold Ratio')
                        ylim([0.5 3])
                        xlim([0 1])
                        subplot(3,3,8)
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_avg =  mean(good_thresh(both_good_ind,:),1);
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_sem =  std(good_thresh(both_good_ind,:),[],1)./sqrt(length(good_ind));
                        errorbar(1:2, mean(good_thresh(both_good_ind,:),1), std(good_thresh(both_good_ind,:),[],1)./sqrt(length(good_ind)),'o')
                        xlim([0.5 2.5])
                        subplot(3,3,9)
                        thresh_ratio = good_thresh(both_good_ind,2)./good_thresh(both_good_ind,1);
                        errorbar(1, mean(thresh_ratio,1), std(thresh_ratio,[],1)./sqrt(length(both_good_ind)),'o')
                        ylim([0.5 3])
                        xlim([0.5 1.5])
                        suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW' ' BaseCon' num2str(pv.basecon_mat(icon))]);
                        fname = ['ThresholdSummary_i' num2str(pv.mouse_mat(imouse)) '_BaseCon' num2str(pv.basecon_mat(icon)) '_' num2str(pv.power_mat(ipow)) 'mW_' pv.task_mat(itask,:) '.pdf'];
                        fn_out = fullfile(rc.fitOutputSummary, fname);
                        exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                    end
                    if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
%                         figure;
%                         x = size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1);
%                         subplot(2,4,5); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1,2),'o'); title('LED'); xlim([x-14 x])
%                         subplot(2,4,1); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1,1),'o'); title('fit'); xlim([x-14 x])
%                         subplot(2,4,2); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2,1),'o'); title('early'); xlim([x-14 x])
%                         subplot(2,4,6); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2,2),'o'); title('LED'); xlim([x-14 x])
%                         subplot(2,4,3); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,3,1),'o'); title('max'); xlim([x-14 x])
%                         subplot(2,4,7); plot(1:x, mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,3,2),'o'); title('LED'); xlim([x-14 x])
%                         subplot(2,4,4); plot(1:x, sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,:,1),2),'o'); title('sum'); xlim([x-14 x])
%                         subplot(2,4,8); plot(1:x, sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,:,2),2),'o'); title('LED'); xlim([x-14 x])
%                         suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
%                         fn_out = fullfile(rc.fitOutputSummary, ['FitEarlyMax_i' num2str(pv.mouse_mat(imouse)) '_BaseCon' num2str(pv.basecon_mat(icon)) '_' num2str(pv.power_mat(ipow)) 'mW_' pv.task_mat(itask,:) '.pdf']]);
%                         exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
%                        figure;
%                         for las = 1:2
%                             fit_early = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,1:2,las),2))==2);
%                             fit_max = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,[1 3],las),2))==2);
%                             early_max = find(squeeze(sum(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max(:,2:3,las),2))==2);
%                             subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1),fit_early,'o'); title('fit_early')
%                             subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), fit_max,'o'); title('fit_max')
%                             subplot(2,3,1+((las-1)*3)); plot(1:size(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).fit_early_max,1), early_max,'o'); title('early_max')
%                         end
%                         suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
%                         fn_out = fullfile(rc.fitOutputSummary, ['FitEarlyMax_pairs_i' num2str(pv.mouse_mat(imouse)) '_'  pv.task_mat(itask,:) '_' num2str(pv.power_mat(ipow)) 'mW.pdf']);
%                         exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                    end
                end
                if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)>0
                    figure;
                    for icon = 1:length(pv.basecon_mat)
                        if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
                            errorbar(pv.basecon_mat(icon), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_avg(:,1), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_sem(:,1), 'ok')
                            hold on
                            errorbar(pv.basecon_mat(icon), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_avg(:,2), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).thresh_sem(:,2), 'oc')
                            hold on
                            ylabel('Threshold (deg)')
                            xlabel('Base contrast')
                            ylim([0 20])
                            suptitle([num2str(pv.mouse_mat(imouse)) ' ' pv.task_mat(itask,:) ' ' num2str(pv.power_mat(ipow)) ' mW']);
                            fname = ['ThresholdSummary_i' num2str(pv.mouse_mat(imouse)) '_' num2str(pv.power_mat(ipow)) 'mW_' pv.task_mat(itask,:) '.pdf'];
                            fn_out = fullfile(rc.fitOutputSummary, fname);
                            exportfig_print(gcf, fn_out, 'FileFormat', 'pdf');
                        end
                    end        
                end
            end
        end
    end
end
x = date;
fn_out = fullfile(rc.fitOutputSummary, ['PVcon_summarystruct_Task' num2str(tasks) '.mat']);
save(fn_out, 'mouse');