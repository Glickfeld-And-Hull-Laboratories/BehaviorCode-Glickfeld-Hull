function plotAllChr2Data(tasks)
mouse = createExptStructChr2(tasks);
pv = behavParamsPV;
rc = behavConstsHADC8;

%plot postion by mouse
figure;
for itask = tasks
    if itask == 1
        subplot(2,2,1)
        xpos = [1 2; 3 4; 5 6];
        pow = 0.25;
        ipow = find(pv.power_mat == pow);
        start = 1;
        for ipos = [1 3 4];
            ns = [];
            for imouse = find(pv.chr2_mat == 1)
                las_mat = [];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                    las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
                                    semilogy(xpos(start,:), [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh], pv.col_str(imouse,:))
                                    hold on
                                end
                            end
                        end
                    end
                end
                if length(las_mat)>0
                    las_avg = mean(las_mat,1);
                    las_sem = std(las_mat,[],1)./sqrt(size(las_mat,1));
                    errorbar(xpos(start,:), las_avg, las_sem, ['o' pv.col_str(imouse,:)]);
                    ns = [ns size(las_mat,1)];
                    mouse(imouse).task(itask).pos(ipos).las_mat= las_mat;
                    mouse(imouse).task(itask).pos(ipos).las_mat_avg = las_avg;
                    mouse(imouse).task(itask).pos(ipos).las_mat_sem = las_sem;
                    mouse(imouse).task(itask).pos(ipos).n = size(las_mat,1);
                else
                    mouse(imouse).task(itask).pos(ipos).las_mat= [];
                    mouse(imouse).task(itask).pos(ipos).las_mat_avg = [];
                    mouse(imouse).task(itask).pos(ipos).las_mat_sem = [];
                    mouse(imouse).task(itask).pos(ipos).n = [];
                end
            end
            text(xpos(start,1), .02, num2str(ns))
            hold on
            start = start+1;
        end
        xlim([0 7])
        ylim([.01 .75])
        ylabel('Threshold Contrast (%)')
        text(1.4,.45, 'Contra')
        text(3.3, .45, 'Middle')
        text(5.2, .45, 'Ipsi')
        title(['Chr2 Task ' num2str(itask) ' by Position All Expt'])

        subplot(2,2,2)
        xpos = [1 2; 3 4; 5 6];
        start = 1;
        for ipos = [1 3 4];
            for imouse = find(pv.chr2_mat == 1)
                las_mat = [];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                    las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
                                end
                            end
                        end
                    end
                end
                if length(las_mat>0)
                    las_mat_norm = las_mat(:,2)./las_mat(:,1);
                    las_mat_avg = mean(las_mat_norm,1);
                    las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat,1));
                    semilogy(xpos(start,1), las_mat_avg, ['o' pv.col_str(imouse,:)]);
                    hold on
                    errorbar(xpos(start,1), las_mat_avg, las_mat_sem, [pv.col_str(imouse,:)]);
                    hold on
                end
            end
            start = start+1;
        end
        ylim([0.5 4])
        xlim([0 7])
        ylabel('Change in Threshold (with PV Activation)')
        text(.8,.55, 'Contra')
        text(2.8, .55, 'Middle')
        text(4.8, .55, 'Ipsi')
        title(['Chr2 Task ' num2str(itask) ' by Position Avg'])

        for ipos = [1 3 4]
            las_mat_all = [];
            n_mice = 0;
            for imouse = find(pv.chr2_mat == 1)
                if mouse(imouse).task(itask).pos(ipos).n > 1
                    las_mat_all = [las_mat_all; mouse(imouse).task(itask).pos(ipos).las_mat];
                    n_mice = n_mice+1;
                end
            end
            las_mat_all_norm = las_mat_all(:,2)./las_mat_all(:,1);
            stats(ipos).las_mat_avg = mean(las_mat_all_norm,1);
            stats(ipos).las_mat_sem = std(las_mat_all_norm,[],1)./sqrt(size(las_mat_all,1));
            stats(ipos).n_sessions = size(las_mat_all,1);
            stats(ipos).n_mice = n_mice;
            las_mat_all_log = log2(las_mat_all);
            stats(ipos).wilkox = signrank(las_mat_all_log(:,1),las_mat_all_log(:,2), 'tail', 'left');
            [H stats(ipos).ks] = kstest2(las_mat_all_log(:,1),las_mat_all_log(:,2), 'Tail', 'larger');
        end
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(tasks) '_mouse.mat']);
        save(pn, 'mouse');
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(tasks) '_stats.mat']);
        save(pn, 'stats');
    end
            
    if itask > 1
    pow = [0.25, NaN, 0.25, NaN, .15];
    subplot(2,2,1)
        xpos = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];
        ipos = 1;
        for icon = 1:length(pv.basecon_mat)
            ns = [];
            y2 = .1;
            for imouse = find(pv.chr2_mat == 1)
                las_mat = [];
                if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                    ipow = find(pv.power_mat == pow(imouse));
                	if ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                            las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh]];
                                            semilogy(xpos(icon,:), [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh], pv.col_str(imouse,:))
                                            hold on
                                        end
                                    end
                                end
                            end
                        end
                        if length(las_mat)>0
                            las_avg = mean(las_mat,1);
                            las_sem = std(las_mat,[],1)./sqrt(size(las_mat,1));
                            errorbar(xpos(icon,:), las_avg, las_sem, ['o' pv.col_str(imouse,:)]);
                            hold on
                            ns = [ns size(las_mat,1)];
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat= las_mat;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat_avg = las_avg;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat_sem = las_sem;
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).n = size(las_mat,1);
                        else
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat= [];
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat_avg = [];
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat_sem = [];
                            mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).n = [];
                        end
                        if itask == 2
                            y = 3;
                            y2 = 2;
                        elseif itask == 3
                            y = 0.1;
                            y2 = 0.05;
                        elseif itask == 4
                            y = 0.1;
                            y2 = 0.05;
                        end
                        text(xpos(icon,1)-.2,y, [num2str(pv.basecon_mat(icon))])
                    end
                end
            end
            text(xpos(icon,1),y2, num2str(ns))
            hold on
        end
        if itask == 2
        xlim([0 11])
        ylim([1 90])
        ylabel('Threshold Orientation (deg)') 
        title(['Chr2 Task ' num2str(itask) ' by Position All Expt'])
        elseif itask == 3
        xlim([0 13])
        ylim([.01 .5])
        ylabel('Threshold Contrast (%)') 
        title(['Chr2 Task ' num2str(itask) ' by Position All Expt'])
        elseif itask == 4
        xlim([0 11])
        ylim([.01 .5])
        ylabel('Threshold Contrast (%)') 
        title(['Chr2 Task ' num2str(itask) ' by Position All Expt'])
        end
        subplot(2,2,2)
        xpos = [1 2; 3 4; 5 6; 7 8; 9 10];
        ipos = 1;
        for imouse = find(pv.chr2_mat == 1)
            if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                ipow = find(pv.power_mat == pow(imouse));
                for icon = 1:length(pv.basecon_mat)
                    las_mat = [];
                    for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                        	if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                        las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh]];
                                        if pv.basecon_mat(icon) == 0
                                            semilogy(log10(.001), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh/mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh, ['s' pv.col_str(imouse,:)]);
                                            hold on
                                        else
                                            semilogy(log10(pv.basecon_mat(icon)), mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh/mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh, ['s' pv.col_str(imouse,:)]);
                                            hold on
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if length(las_mat>0)
                        las_mat_norm = las_mat(:,2)./las_mat(:,1);
                        las_mat_avg = mean(las_mat_norm,1);
                        las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat,1));       
                        if pv.basecon_mat(icon) == 0
                            semilogy(log10(.001), las_mat_avg, pv.col_str(imouse,:));
                            hold on
                            errorbar(log10(.001), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                        else
                            semilogy(log10(pv.basecon_mat(icon)), las_mat_avg, pv.col_str(imouse,:));
                            hold on
                            errorbar(log10(pv.basecon_mat(icon)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                        end
                        hold on
                    end
                end
            end
        end
        ylim([.5 5])
        xlim([-3 .2])
        ylabel('Change in Threshold (with PV Activation)')
        xlabel('Base contrast')
        title(['Chr2 Task ' num2str(itask) ' by Base Contrast Avg'])
        
        for icon = 1:length(pv.basecon_mat)
            ipos = 1;
            las_mat_all = [];
            n_mice = 0;
            for imouse = find(pv.chr2_mat == 1)
                pow = [0.25, NaN, 0.25, NaN, .15];
                ipow = find(pv.power_mat == pow(imouse));
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind > 0
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).n > 0
                        las_mat_all = [las_mat_all; mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).las_mat];
                        n_mice = n_mice+1;
                    end
                end
            end
            if size(las_mat_all,1)>1
                las_mat_all_norm = las_mat_all(:,2)./las_mat_all(:,1);
                stats(icon).las_mat_avg = mean(las_mat_all_norm,1);
                stats(icon).las_mat_sem = std(las_mat_all_norm,[],1)./sqrt(size(las_mat_all,1));
                stats(icon).n_sessions = size(las_mat_all,1);
                stats(icon).n_mice = n_mice;
                las_mat_all_log = log2(las_mat_all);
                stats(icon).wilkox = signrank(las_mat_all_log(:,1),las_mat_all_log(:,2), 'tail', 'left');
                [H stats(icon).ks] = kstest2(las_mat_all_log(:,1),las_mat_all_log(:,2), 'Tail', 'larger');
            else
                stats(icon).las_mat_avg = [];
                stats(icon).las_mat_sem = [];
                stats(icon).n_sessions = [];
                stats(icon).n_mice = n_mice;
                stats(icon).wilkox = [];
                stats(icon).ks = [];
            end
        end
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(tasks) '_stats.mat']);
        save(pn, 'stats');
    end    
    if itask == 3
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-summary-' date '.pdf']);
        exportfig_print(gcf, pn, 'FileFormat', 'pdf');
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(tasks) '_mouse.mat']);
        save(pn, 'mouse');
    end
    if itask == 4;
        subplot(2,2,3)
        for imouse = find(pv.chr2_mat == 1)
            ipos = 1;
            if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                pow = [.25, NaN, 0.25, NaN, .15];
                ipow = find(pv.power_mat == pow(imouse));
                for icon = 1:length(pv.basecon_mat)
                    thresh_mat = [];
                    if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                        thresh = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh;
                                        thresh_mat = [thresh_mat; thresh];
                                        plot(log(pv.basecon_mat(icon)), thresh, 'ob')
                                        hold on
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        xlabel('Base Contrast')
        ylabel('Contrast Threshold (%)')
        xlim([-3 0.2])
        ylim([.1 .5])
        subplot(2,2,4)
        con_col = strvcat('k', 'b', 'c', 'g');
        for imouse = find(pv.chr2_mat == 1)
            ipos = 1;
            if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                pow = [.25, NaN, 0.25, NaN, .15];
                ipow = find(pv.power_mat == pow(imouse));
                for icon = 1:length(pv.basecon_mat)
                    thresh_mat = [];
                    if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                            thresh = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh;
                                            thresh_ratio = mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh;
                                            plot(thresh, thresh_ratio, ['o' con_col(icon,:)])
                                            hold on
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        xlabel('Threshold')
        ylabel('Change in Threshold')
        xlim([0 .5])
        ylim([.5 3])
        pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-summary-' date '.pdf']);
        exportfig_print(gcf, pn, 'FileFormat', 'pdf');        
    end
%avg by position by mouse


% %by power by mouse
% subplot(2,2,3)
% start = [1 2];
% xpos = [];
% for rep = 1:length(pv.power_mat)
%     xpos = [xpos; start + (rep-1)*[2 2]];
% end
% ipos = 1;
% for imouse = find(pv.chr2_mat == 1)
%     for ipow = 1:length(pv.power_mat)
%         las_mat = [];
%         for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
%             if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh>0
%                 las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
%                 semilogy(xpos(ipow,:), [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh], ...
%                     pv.col_str(imouse,:));
%                 hold on
%             end
%         end
%         if length(las_mat>0)
%             las_mat_avg = mean(las_mat,1);
%             las_mat_sem = std(las_mat,[],1)./sqrt(size(las_mat,1));
%             errorbar(xpos(ipow,:), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
%             hold on
%             text(xpos(ipow,1), 3, num2str(pv.power_mat(ipow)));
%         end
%     end
% end
% if itask == 1
%     ylim([0.01 0.5])
%     ylabel('Threshold contrast (%)')
% else
%     ylim([1 50])
%     ylabel('Threshold orientation (deg)')
% end
% title(['Chr2 Task ' num2str(itask) ' by Power All Expt'])

%change in thresh by power by mouse
area_mm2 = [0.7802; NaN; 0.5453; NaN; 0.4016];
if itask < 3
    for imouse = find(pv.chr2_mat == 1)
        las_mat_all = [];
        for ipow = 1:length(pv.power_mat)
            ipos = 1;
            las_mat = [];
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                	if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
                            end
                        end
                    end
                end
            end
            if length(las_mat>0)
                las_mat_norm = las_mat(:,2)./las_mat(:,1);
                las_mat_avg = mean(las_mat_norm,1);
                las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                subplot(2,2,3)
                semilogy(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, ['o' pv.col_str(imouse,:)]);
                hold on
                errorbar(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, las_mat_sem, [ pv.col_str(imouse,:)]);
                hold on
                mouse(imouse).task(itask).pos(ipos).pow(ipow).las_mat_avg = las_mat_avg;
                mouse(imouse).task(itask).pos(ipos).pow(ipow).las_mat_sem = las_mat_sem;
            end
        end
    end
    subplot(2,2,3)
    ylim([0.5 4])
    xlim([-3 1])
    xlabel('log(laser power)')
    ylabel('Change in Threshold')
    title(['Chr2 Task ' num2str(itask) ' by Power Avg'])
    %con = 1
    if itask > 1
        for ipow = 1:length(pv.power_mat)
            ns = [];
            for imouse = find(pv.chr2_mat == 1)
                ipos = 1;
                icon = 1;
                las_mat = [];
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                    las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh]];
                                    semilogy(log10(pv.power_mat(ipow)./area_mm2(imouse)), (mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh./mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh), ['s' pv.col_str(imouse,:)]);
                                    hold on
                                end
                            end
                        end
                    end
                end
                if length(las_mat>0)
                    las_mat_norm = las_mat(:,2)./las_mat(:,1);
                    las_mat_avg = mean(las_mat_norm,1);
                    las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));   
                    subplot(2,2,4)
                    semilogy(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, ['o' pv.col_str(imouse,:)]);
                    hold on
                    errorbar(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, las_mat_sem, pv.col_str(imouse,:));
                    length(las_mat)
                    ns = [ns size(las_mat,1)];
                end
            end
            text(log10(pv.power_mat(ipow)./area_mm2(imouse)),0.3,num2str(ns));
            hold on
        end
        subplot(2,2,4)
        ylim([.5 5])
        xlim([-3 1])
        xlabel('laser power')
        ylabel('Change in Threshold')
        title(['Chr2 Task ' num2str(itask) ' by Power when Con = 1'])
    end
    pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(tasks) '_mouse.mat']);
    save(pn, 'mouse');
    pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-summary-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
    %max estimate
    figure;
    for imouse = find(pv.chr2_mat == 1)
        for ipow = 1:length(pv.power_mat)
            ipos = 1;
            las_mat = [];
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early == 1
                	if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).miss == 1
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).max mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).max]];
                            end
                        end
                    end
                end
            end
            if length(las_mat>0)
                las_mat_norm = las_mat(:,2);
                las_mat_avg = mean(las_mat_norm,1);
                las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                subplot(2,2,1)
                errorbar(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                hold on
                las_mat_norm = las_mat(:,2)./las_mat(:,1);
                las_mat_avg = mean(las_mat_norm,1);
                las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                subplot(2,2,2)
                errorbar(log10(pv.power_mat(ipow)./area_mm2(imouse)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                hold on
            end
        end
    end
    subplot(2,2,1)
    ylim([.75 1.25])
    xlim([-3 0])
    xlabel('log(laser power)')
    ylabel('Max with laser')
    subplot(2,2,2)
    ylim([.75 1.25])
    xlim([-3 0])
    xlabel('log(laser power)')
    ylabel('Change in Max')
    suptitle(['Chr2 Task ' num2str(itask) ' by Power Avg'])

    if itask == 2
    pow = [.25, NaN, 0.25, NaN, .15];
        xpos = [1 2; 3 4; 5 6];
        ipos = 1;
        for imouse = find(pv.chr2_mat == 1)
            if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                ipow = find(pv.power_mat == pow(imouse));
                for icon = 1:length(pv.basecon_mat)
                    las_mat = [];
                    if ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                        las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).max mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).max]];
                                        hold on
                                        end
                                    end
                                end
                            end
                        end
                        if length(las_mat)>0
                            las_avg = mean(las_mat(:,2),1);
                            las_sem = std(las_mat(:,2),[],1)./sqrt(size(las_mat(:,2),1));
                            subplot(2,2,3)
                            errorbar(log10(pv.basecon_mat(icon)), las_avg, las_sem, ['o' pv.col_str(imouse,:)]);
                            hold on
                            las_mat_norm = las_mat(:,2)./las_mat(:,1);
                            las_avg = mean(las_mat_norm,1);
                            las_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                            subplot(2,2,4)
                            errorbar(log10(pv.basecon_mat(icon)), las_avg, las_sem, ['o' pv.col_str(imouse,:)]);
                            hold on
                        end
                    end
                end
            end
        end
        subplot(2,2,3)
        xlim([-1 0])
        ylim([.75 1.25])
        ylabel('Max with laser') 
        xlabel('log(Contrast)')
        subplot(2,2,4)
        xlim([-1 0])
        ylim([.75 1.25])
        ylabel('Change in Max') 
        xlabel('log(Contrast)')
    end
    pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-max-saturation-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

    %for task 2 only
    if itask == 2
        las_con_mat = zeros(2,length(pv.mouse_mat), length(pv.basecon_mat));
        for imouse = find(pv.chr2_mat == 1)
            ipos = 1;
            if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                figure;
                pow = [.25, NaN, 0.25, NaN, .15];
                ipow = find(pv.power_mat == pow(imouse));
                for icon = 1:length(pv.basecon_mat)
                    las_mat = [];
                    if ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).early == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).miss == 1
                                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).good == 1
                                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).good == 1
                                            las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).expt(iexp).laser(2).thresh]];
                                        end
                                    end
                                end
                            end
                        end
                        if length(las_mat)>0
                            las_avg = mean(las_mat,1);
                            las_sem = std(las_mat,[],1)./sqrt(size(las_mat,1));
                            las_con_mat(:,imouse, icon) = las_avg;
                            subplot(2,2,1)
                            plot(log10(pv.basecon_mat(icon))*ones(size(las_mat,1),1), (las_mat(:,1)), 'ok')
                            hold on
                            plot(log10(pv.basecon_mat(icon))*ones(size(las_mat,1),1), (las_mat(:,2)), 'or')
                            xlim([-1 .25])
                            ylim([0 80])
                            subplot(2,2,2)
                            errorbar((pv.basecon_mat(icon)), (las_avg(:,1)), (las_sem(:,1)), 'ok');
                            hold on
                            errorbar((pv.basecon_mat(icon)), (las_avg(:,2)), (las_sem(:,2)), 'or');
                            xlim([0 1.1])
                            ylim([0 80])
                            subplot(2,2,3)
                            errorbar(log10(pv.basecon_mat(icon)), (las_avg(:,1)), (las_sem(:,1)), 'ok');
                            hold on
                            errorbar(log10(pv.basecon_mat(icon)), (las_avg(:,2)), (las_sem(:,2)), 'or');
                            xlim([-1 .25])
                            ylim([0 80])
                            ylabel('Threshold Orientation (deg)')
                            xlabel('Contrast')
                            subplot(2,2,4)
                            plot([pv.basecon_mat(icon) pv.basecon_mat(icon)], [las_avg(:,1)-las_sem(:,1) las_avg(:,1)+las_sem(:,1)], '-k')
                            hold on
                            plot(pv.basecon_mat(icon), las_avg(:,1), 'ok');
                            hold on
                            plot([pv.basecon_mat(icon) pv.basecon_mat(icon)], [las_avg(:,2)-las_sem(:,2) las_avg(:,2)+las_sem(:,2)], '-r')
                            hold on
                            plot(pv.basecon_mat(icon), las_avg(:,2), 'or');
                            hold on
                            set(gca, 'xscale', 'log');
                            set(gca, 'yscale', 'log');
                            xlim([.1 1.1])
                            ylim([1 80])
                        end
                    end
                end
                suptitle(['Chr2 task ' num2str(itask) ' mouse i' num2str(pv.mouse_mat(imouse))]) 
                pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-threshold-summary-i' num2str(pv.mouse_mat(imouse)) '-' date '.pdf']);
                exportfig_print(gcf, pn, 'FileFormat', 'pdf');
            end
        end
    end
    if tasks == 1:2;
        if itask == 2;
            figure;
            for imouse = find(pv.chr2_mat == 1)
                ipos = 1;
                if length(mouse(imouse).task(itask).pos(ipos).ind)>0
                    pow = [.25, NaN, 0.25, NaN, .15];
                    ipow = find(pv.power_mat == pow(imouse));
                    apparent_con = 1-(1./mouse(imouse).task(itask).pos(ipos).pow(ipow).las_mat_avg);
                    for icon = 1:length(pv.basecon_mat)
                        if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).con(icon).ind)>0
                            if size(las_con_mat(:,imouse,icon),1)>0
                                if length(pv.basecon_mat)>= icon+2
                                    step = icon+1;
                                    con_change = (pv.basecon_mat(icon+1))./pv.basecon_mat(icon);
                                else
                                    step = icon-1;
                                    con_change = (pv.basecon_mat(icon))./pv.basecon_mat(icon-1);
                                end
                                act_ori = squeeze(las_con_mat(2,imouse, icon));
                                ori_diff = abs(squeeze(las_con_mat(1,imouse, step))-squeeze(las_con_mat(1, imouse, icon)));
                                pred_ori = squeeze(las_con_mat(1, imouse, icon)) + (apparent_con./con_change)*ori_diff;
                                plot(act_ori, pred_ori, [pv.mark_str(icon,:) pv.col_str(imouse,:)]);
                                hold on
                            end
                        end
                    end
                end
            end
            x = 1:1:90;
            y = x;
            plot(x,y,'--k')
            xlim([0 90])
            ylim([0 90])
            xlabel('Actual orientation threshold (deg)')
            ylabel('Predicted orientation threshold (deg)')
            suptitle(['Chr2 task ' num2str(itask) ' predicted ori shift']) 
            pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-predicted-ori-shift-' date '.pdf']);
            exportfig_print(gcf, pn, 'FileFormat', 'pdf');
        end
    end
end
end
end
