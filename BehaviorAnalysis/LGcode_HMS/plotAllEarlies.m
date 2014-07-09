function plotAllEarlies(itask)

mouse = createExptStructPV(itask);
pv = behavParamsPV;
rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);


for imouse = 1:length(pv.mouse_mat)
    for ipos = 1:length(pv.pos_mat)
        for ipow = 1:length(pv.power_mat)
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1;
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1;
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1;
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1;
                                fn = fullfile('/Users/lindsey/Desktop/Data', ['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp))) '.mat']);
                                load(fn);
                                trials = cell2mat(xd.TrialRangeToUse(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp)));
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
                                block1_ind_early = find(RTs(block1_ind_good)<=100);
                                block2_ind_early = find(RTs(block2_ind_good)<=100);
                                block1_ind_miss = find(RTs(block1_ind_good)>550);
                                block2_ind_miss = find(RTs(block2_ind_good)>550);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all = length(block1_ind_good);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all = length(block2_ind_good);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early = length(block1_ind_early);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early = length(block2_ind_early);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_miss = length(block1_ind_miss);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_miss = length(block2_ind_miss);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_correct = length(block1_ind_good)-length(block1_ind_early)-length(block1_ind_miss);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_correct = length(block2_ind_good)-length(block2_ind_early)-length(block2_ind_miss);
                                block1_early_prevRT = RTs(block1_ind_early(2:end)-1);
                                block2_early_prevRT = RTs(block2_ind_early(2:end)-1);
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_prev_early = length(find(block1_early_prevRT<125)); 
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_prev_early = length(find(block2_early_prevRT<125));
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_prev_miss = length(find(block1_early_prevRT>550)); 
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_prev_miss = length(find(block2_early_prevRT>550));
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_prev_correct = length(block1_ind_early)-length(find(block1_early_prevRT>550))-length(find(block1_early_prevRT<125)); 
                                mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_prev_correct = length(block2_ind_early)-length(find(block2_early_prevRT>550))-length(find(block2_early_prevRT<125));
                            end
                            if itask > 1
                                for icon = 1:length(pv.basecon_mat)
                                    for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                                        fn = fullfile('/Users/lindsey/Desktop/Data', ['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp))) '.mat']);
                                        load(fn);
                                        trials = cell2mat(xd.TrialRangeToUse(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp)));
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
                                        block1_ind_early = find(RTs(block1_ind_good)<=100);
                                        block2_ind_early = find(RTs(block2_ind_good)<=100);
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_all = length(block1_ind_good);
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_all = length(block2_ind_good);
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_early = length(block1_ind_early);
                                        mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_early = length(block2_ind_early);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%avg change in earlies by power by mouse
figure;
subplot(2,2,1)
for imouse = 1:length(pv.mouse_mat)
    for ipos = 1
        for ipow = 1:length(pv.power_mat)
            las_mat = [];
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                                pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all;
                                pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all;
                                las_mat = [las_mat; pct_earlies_block1 pct_earlies_block2];
                            end
                        end
                    end
                end
            end
            if length(las_mat>0)
                las_mat_norm = las_mat(:,2)./las_mat(:,1);
                las_mat_avg = mean(las_mat_norm,1);
                las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                errorbar(log10(pv.power_mat(ipow)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                hold on
            end
        end
    end
end
ylim([0.5 1.5])
xlim([floor(log10(pv.power_mat(end))) ceil(log10(pv.power_mat(1)))])
xlabel('log(laser power)')
ylabel('Change in Earlies')
title(['Task ' num2str(itask) ' Earlies by power'])

% change in earlies by power by mouse
subplot(2,2,2)
for imouse = 1:length(pv.mouse_mat)
    for ipos = 1
        for ipow = 1:length(pv.power_mat)
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                                pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all;
                                pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all;
                                plot(log10(pv.power_mat(ipow)), pct_earlies_block2./pct_earlies_block1, ['o' pv.col_str(imouse,:)])
                                hold on
                            end
                        end
                    end
                end
            end
        end
    end
end
ylim([0.5 1.5])
xlim([floor(log10(pv.power_mat(end))) ceil(log10(pv.power_mat(1)))])
xlabel('log(laser power)')
ylabel('Change in Earlies')
title(['Task ' num2str(itask) ' Earlies by power'])

% change in earlies by change in thresh by mouse
subplot(2,2,3)
for imouse = 1:length(pv.mouse_mat)
    for ipos = 1
        for ipow = 1:length(pv.power_mat)
            las_mat = [];
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                                pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all;
                                pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all;
                                change_thresh =  mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh;
                                las_mat = [las_mat; change_thresh pct_earlies_block2./pct_earlies_block1];
                            end
                        end
                    end
                end
            end
            if size(las_mat,1)>0
                las_mat_avg = mean(las_mat,1);
                if size(las_mat,1)>1
                    las_mat_sem = std(las_mat,1)./sqrt(size(las_mat,1));
                else
                    las_mat_sem = [0 0];
                end
                errorbarxy(las_mat_avg(1), las_mat_avg(2), las_mat_sem(1), las_mat_sem(2),[],[],['o' pv.col_str(imouse,:)], pv.col_str(imouse,:));
                hold on
            end
        end
    end
end
ylim([.5 1.5])
xlim([.5 3.5])
xlabel('Change in Thresh')
ylabel('Change in Earlies')
title(['Task ' num2str(itask) ' Earlies by Change in Threshold'])
pn = fullfile(rc.fitOutputSummary, ['Task' num2str(itask) '-earlies-summary-' date '.pdf']);
exportfig_print(gcf, pn, 'FileFormat', 'pdf');

% %earlies by contrast
% if itask > 1
%     subplot(2,2,4)
%     for imouse = 1:length(pv.mouse_mat)
%         for ipos = 1
%             for icon = 1:length(pv.basecon_mat)
%                 for ipow = 1:length(pv.power_mat)
%                     las_mat = [];
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
%                             if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
%                                 pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_all;
%                                 pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_all;
%                                 las_mat = [las_mat; pct_earlies_block1 pct_earlies_block2];
%                             end
%                         end
%                     end
%                     if length(las_mat>0)
%                         las_mat_norm = las_mat(:,2)./las_mat(:,1);
%                         las_mat_avg = mean(las_mat_norm,1);
%                         las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
%                         errorbar(log10(pv.basecon_mat(icon)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
%                         hold on
%                     end
%                 end
%             end
%         end
%     end
%     ylim([0.5 1.5])
%     xlim([floor(log10(pv.basecon_mat(end))) ceil(log10(pv.basecon_mat(1)))])
%     xlabel('log(base contrast)')
%     ylabel('Change in Earlies')
%     title(['Task ' num2str(itask) ' Earlies by base contrast'])
%     pn = fullfile(rc.fitOutputSummary, ['Task' num2str(itask) '-earlies-summary-' date '.pdf']);
%     exportfig_print(gcf, pn, 'FileFormat', 'pdf');    
% 
% %block1 earlies by contrast
% figure;
%     for imouse = 1:length(pv.mouse_mat)
%         for ipos = 1
%             for icon = 1:length(pv.basecon_mat)
%                 for ipow = 7
%                     early_mat = [];
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                             pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_all;
%                             early_mat = [early_mat; pct_earlies_block1];
%                         end
%                     end
%                     if length(early_mat>0)
%                         early_mat_avg = mean(early_mat,1);
%                         early_mat_sem = std(early_mat,[],1)./sqrt(size(early_mat,1));
%                         errorbar(log10(pv.basecon_mat(icon)), early_mat_avg, early_mat_sem, ['o' pv.col_str(imouse,:)]);
%                         hold on
%                     end
%                 end
%             end
%         end
%     end
% %     xlim([floor(log10(pv.basecon_mat(end))) ceil(log10(pv.basecon_mat(1)))])
%     xlabel('log(base contrast)')
%     ylabel('Block1 Early %')
%     title(['Task ' num2str(itask) ' Block 1 earlies by base contrast'])
% end
% 
% %timecourse of earlies/threshold
%     for imouse = 1:length(pv.mouse_mat)
%         ipos = 1;
%         if length(mouse(imouse).task(itask).pos(ipos).ind)>0
%             for ipow = 1:length(pv.power_mat)
%                 pow_exist(ipow) = ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind);
%             end
%             max_pow = pv.power_mat(min(find(pow_exist == 1)));
%             figure;
%             if itask>1    
%                 for ipow = 1:length(pv.power_mat)
%                     for icon = 1:length(pv.basecon_mat)
%                         con = pv.basecon_mat(icon);
%                         subplot(2,2,1)
%                         cool_map = cool(64);
%                         col = ceil((pv.power_mat(ipow)./max_pow)*64);
%                         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                             if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
%                                 semilogy(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh,'o', 'color', cool_map(col,:))
%                                 hold on
%                             end
%                         end
%                         subplot(2,2,2)
%                         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                             if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
%                                 if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
%                                     plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp), (mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh)./(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh), 'o', 'color', cool_map(col,:))
%                                     hold on
%                                 end
%                             end
%                         end
%                         subplot(2,2,3)
%                         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                             if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
%                                 if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
%                                     pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_all;
%                                     pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_all;
%                                     delta_earlies = pct_earlies_block2./pct_earlies_block1;
%                                     plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp), delta_earlies, 'o', 'color', cool_map(col,:))
%                                     hold on
%                                 end
%                             end
%                         end
%                         subplot(2,2,4)
%                         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
%                             if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
%                                 if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
%                                     total_earlies = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_early+mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_early;
%                                     total_trials = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).n_all+mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).n_all;
%                                     early_rate = total_earlies./total_trials;
%                                     plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind(iexp), early_rate, 'o', 'color', cool_map(col,:))
%                                     hold on
%                                 end
%                             end
%                         end
%                     end
%                 end
%             elseif itask ==1
%                 for ipow = 1:length(pv.power_mat)
%                     subplot(2,2,1)
%                     cool_map = cool(64);
%                     col = ceil((pv.power_mat(ipow)./max_pow)*64);
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                             semilogy(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp), mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh,'o', 'color', cool_map(col,:))
%                             hold on
%                         end
%                     end
%                     subplot(2,2,2)
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                         	if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
%                                 plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp), (mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh)./(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh), 'o', 'color', cool_map(col,:))
%                                 hold on
%                             end
%                         end
%                     end
%                     subplot(2,2,3)
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                             pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all;
%                             pct_earlies_block2 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all;
%                             delta_earlies = pct_earlies_block2./pct_earlies_block1;
%                             plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp), delta_earlies, 'o', 'color', cool_map(col,:))
%                             hold on
%                         end
%                     end
%                     subplot(2,2,4)
%                     for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%                         if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
%                             total_earlies = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early+mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early;
%                             total_trials = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all+mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all;
%                             early_rate = total_earlies./total_trials;
%                             plot(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp), early_rate, 'o', 'color', cool_map(col,:))
%                             hold on
%                         end
%                     end
%                 end
%             end
%             subplot(2,2,1)
%             xlabel('Day')
%             ylabel('Baseline Threshold')
%             ylim([.01 .5])
%             if itask==2
%                 ylim([1 50])
%             end
%             title('Baseline threshold')
%             subplot(2,2,2)
%             xlabel('Day')
%             ylabel('Change in Threshold')
%             ylim([.5 3])
%             title('Change in threshold')
%             subplot(2,2,3)
%             xlabel('Day')
%             ylabel('Change in FAs')
%             ylim([.5 2])
%             title('Change in false alarms')
%             subplot(2,2,4)
%             xlabel('Day')
%             ylabel('% FA')
%             ylim([0 1])
%             title('False alarm rate')           
%             suptitle(['Timecourses for i' num2str(pv.mouse_mat(imouse))])
%             pn = fullfile(rc.fitOutputSummary, ['Task' num2str(itask) '-earlies-summary-i' num2str(pv.mouse_mat(imouse)) '-' date '.pdf']);
%             exportfig_print(gcf, pn, 'FileFormat', 'pdf');  
%         end
%     end

    
%trials preceding false alarms
% figure;
% for imouse = 1
%     for ipos = 1
%         for ipow = 1:length(pv.power_mat)
%             for iblock = 1:2
%             	mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_miss = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_correct = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_early = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_miss = 0;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_correct = 0;
%                 for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%                     if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh > 0
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_miss = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_miss + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_miss];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_correct = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_correct + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_correct];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_early = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_early + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_prev_early];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_miss = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_miss + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_prev_miss];
%                         mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_correct = [mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_correct + mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_prev_correct];
%                     end
%                 end
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_early = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_miss = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_miss./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_correct = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_correct./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_all;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_prev_early = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_prev_miss = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_miss./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early;
%                 mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_prev_correct = mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_prev_correct./mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).n_early;                
%             end
%             if mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(iblock).pct_early > 0
% %                 subplot(2,3,1)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_early, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_early, 'or')
% %                 ylim([0 1])
% %                 title('Early (control)')
% %                 subplot(2,3,2)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_miss, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_miss, 'or')
% %                 ylim([0 1])
% %                 title('Miss (control)')
% %                 subplot(2,3,3)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_correct, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_correct, 'or')
% %                 ylim([0 1])
% %                 title('Correct (control)')
% %                 subplot(2,3,4)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_early, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_early, 'or')
% %                 ylim([0 1])
% %                 title('Early (Chr2)')
% %                 subplot(2,3,5)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_miss, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_miss, 'or')
% %                 ylim([0 1])
% %                 title('Miss (Chr2)')
% %                 subplot(2,3,6)
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_correct, 'ok') 
% %                 hold on
% %                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_correct, 'or')
% %                 ylim([0 1])
% %                 title('Correct (Chr2)')
%                 subplot(2,3,1)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_early, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_early, 'or')
%                 ylim([0 1])
%                 title('Early (control)')
%                 subplot(2,3,2)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_miss, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_miss, 'or')
%                 ylim([0 1])
%                 title('Miss (control)')
%                 subplot(2,3,3)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_correct, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_correct, 'or')
%                 ylim([0 1])
%                 title('Correct (control)')
%                 subplot(2,3,4)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_early, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_early, 'or')
%                 ylim([0 1])
%                 title('Early (Prev)')
%                 subplot(2,3,5)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_miss, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_miss, 'or')
%                 ylim([0 1])
%                 title('Miss (Prev)')
%                 subplot(2,3,6)
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(1).pct_prev_correct, 'ok') 
%                 hold on
%                 plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).laser(2).pct_prev_correct, 'or')
%                 ylim([0 1])
%                 title('Correct (Prev)')
%             end
%         end
%     end
% end
                
end
