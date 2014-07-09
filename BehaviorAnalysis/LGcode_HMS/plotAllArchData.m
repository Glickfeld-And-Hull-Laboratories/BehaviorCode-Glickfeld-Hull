function plotAllArchData(itask)
mouse = createExptStructArch(itask);
pv = behavParamsPV;
rc = behavConstsHADC8;

figure;
%actual threshold by power by mouse (no laser)
subplot(2,2,1)
start = [1 2];
ipos = 1;
for imouse = find(pv.arch_mat == 1)
    for ipow = 1:length(pv.power_mat)
        if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)>0
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).AC)
                    if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh>0
                        semilogy(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh, ['o' pv.col_str(imouse,:)]);
                        hold on
                    end
                end
            end  
        end
    end
end
if itask == 1
    ylim([0.1 1])
    ylabel('Threshold contrast (%)')
else
    ylim([0.1 10])
    ylabel('Threshold orientation (deg)')
end
title(['Baseline threshold by Power All Expt'])

%change in thresh by power by mouse
subplot(2,2,2)
for imouse = find(pv.arch_mat == 1)
    for ipow = 1:length(pv.power_mat)
        las_mat = [];
        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
            if isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).AC)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh>0
                    plot(log10(pv.power_mat(ipow)), mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh, ['o' pv.col_str(imouse,:)]);
                    hold on
                end
            end
        end
    end
end
ylim([0 2])
xlim([floor(log10(pv.power_mat(end))) ceil(log10(pv.power_mat(1)))])
xlabel('log(laser power)')
ylabel('Change in Threshold')
title(['By Power All Expt'])

%avg change in thresh by power by mouse
subplot(2,2,3)
for imouse = find(pv.arch_mat == 1)
    for ipow = 1:length(pv.power_mat)
        las_mat = [];
        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
            if isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).AC)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh>0
                    las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
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
ylim([0 2])
xlim([floor(log10(pv.power_mat(end))) ceil(log10(pv.power_mat(1)))])
xlabel('log(laser power)')
ylabel('Change in Threshold')
title(['By Power Avg'])

suptitle(['Arch Task ' num2str(itask)])
pn = fullfile(rc.fitOutputSummary, ['Arch-Task' num2str(itask) '-DClaser-summary-' date '.pdf']);
exportfig_print(gcf, pn, 'FileFormat', 'pdf');

start = 1;
colstr = strvcat('k', 'c');
for imouse = find(pv.arch_mat == 1)
    for ipos = 1:length(pv.pos_mat)
        if mouse(imouse).task(itask).pos(ipos).ind>0
            pos = pv.pos_mat(ipos);
            figure;
            start = 1;
            for ipow = 1:length(pv.power_mat)
                if mouse(imouse).task(itask).pos(ipos).pow(ipow).b2probe == 1
                    if start > 4
                        fignum = 1;
                        suptitle(['Mouse = i' num2str(pv.mouse_mat(imouse)) ' Pos = ' num2str(pos)]);
                        pn = fullfile(rc.fitOutputSummary, ['Arch-Task' num2str(itask) '-b2probe-summary-mouse-i' num2str(pv.mouse_mat(imouse)) '-pos' num2str(pos) '-fignum' num2str(fignum) '-' date '.pdf']);
                        exportfig_print(gcf, pn, 'FileFormat', 'pdf');
                        figure;
                        start = 1;
                    end
                    subplot(2,2,start)
                    for iblock = 1:2
                        pct_corr = (mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nCorrect)./(mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nMissed+mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nCorrect);
                        plot(1:2, pct_corr, ['o' colstr(iblock,:)]);
                        hold on
                        text(1,.15+(iblock*.1),num2str(mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nMissed(1)+mouse(imouse).task(itask).pos(ipos).pow(ipow).block(iblock).nCorrect(1)), 'Color', colstr(iblock,:))
                    end
                    start = start+1;
                    title(['Power = ' num2str(pv.power_mat(ipow))])
                    ylim([0 1])
                    xlim([0 3])
                end
            end
        end
        suptitle(['Mouse = i' num2str(pv.mouse_mat(imouse)) ' Pos = ' num2str(pos)]);
        pn = fullfile(rc.fitOutputSummary, ['Arch-Task' num2str(itask) '-b2probe-summary-mouse-i' num2str(pv.mouse_mat(imouse)) '-pos' num2str(pos) '-' date '.pdf']);
        exportfig_print(gcf, pn, 'FileFormat', 'pdf');
    end
end
    
    % pulsed laser analysis
    figure;
    xpos = [1 2; 3 4; 5 6];
    for imouse = find(pv.arch_mat == 1)
        las_mat_all = [];
        start = 1;
        for ipow = 1:length(pv.power_mat)
            pass  = 0;
            for ipos = 1:length(pv.pos_mat);
                if length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind>0)
                    las_mat = [];
                    for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                        if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).AC == 1
                            if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).good == 1
                                if mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).good == 1
                                    subplot(2,2,1)
                                    las_mat = [las_mat; [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh]];
                                    semilogy(xpos(start,:), [mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh], pv.col_str(imouse,:));
                                    hold on
                                    pass = 1;
                                end
                            end
                        end
                    end
                    if length(las_mat>0)
                       las_mat_avg = mean(las_mat,1);
                       las_mat_sem = std(las_mat,[],1)./sqrt(size(las_mat,1));
                       errorbar(xpos(start,:), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                       hold on
                       text(xpos(start,1), 0.02, num2str(pv.power_mat(ipow)));
                       hold on
                    end
                    if length(las_mat>0)
                        subplot(2,2,2)
                        las_mat_norm = las_mat(:,2)./las_mat(:,1);
                        las_mat_avg = mean(las_mat_norm,1);
                        las_mat_sem = std(las_mat_norm,[],1)./sqrt(size(las_mat_norm,1));
                        errorbar(log10(pv.power_mat(ipow)), las_mat_avg, las_mat_sem, ['o' pv.col_str(imouse,:)]);
                        hold on
                    end
                end
            end
            if pass == 1;
                start = start+1;
            end
        end
    end
    subplot(2,2,1)
    ylim([0.01 1])
    xlim([0 6])
    ylabel('Contrast threshold')
    title(['Arch Task ' num2str(itask) ' by Power'])
    subplot(2,2,2)
    ylim([0.5 1.5])
    xlim([0 2])
    ylabel('Change in threshold')
    title(['Arch Task ' num2str(itask) ' by Power Avg'])
    suptitle(['Arch Task ' num2str(itask)])
    pn = fullfile(rc.fitOutputSummary, ['Arch-Task' num2str(itask) '-AClaser-summary-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end

