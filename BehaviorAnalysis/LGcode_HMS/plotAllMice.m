function plotAllMice

mouse = createExptStructAll;
pv = behavParamsAll;
rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);


% avg threshold/earlies by mouse

figure;
subplot(1,2,1)
for imouse = 1:length(pv.mouse_mat)
    z =0;
    x_mat = [];
    for itask = [1 3]
        ipos = 1;
        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
            if ~isnan(mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh)
                z = 1;
                x_mat = [x_mat; mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh];
            end
        end
    end
    if z == 1;
    plot(imouse*ones(size(x_mat,1),1), x_mat, 'ok');
    x_mat_avg = mean(x_mat,1);
    x_mat_sem = std(x_mat,1)./sqrt(size(x_mat,1));
    hold on
    errorbar(imouse, x_mat_avg, x_mat_sem, 'or')
    end
end
ylim([0 .3])
ylabel('Threshold')

subplot(1,2,2)
for imouse = 1:length(pv.mouse_mat)
    z =0;
    x_mat = [];
    for itask = [1 3]
        ipos = 1;
        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
            if ~isnan(mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh)
            	pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).expt(iexp).n_early./mouse(imouse).task(itask).pos(ipos).expt(iexp).n_all;
                z = 1;
                x_mat = [x_mat; pct_earlies_block1];
            end
        end
    end
    if z == 1;
    plot(imouse*ones(size(x_mat,1),1), x_mat, 'ok');
    x_mat_avg = mean(x_mat,1);
    x_mat_sem = std(x_mat,1)./sqrt(size(x_mat,1));
    hold on
    errorbar(imouse, x_mat_avg, x_mat_sem, 'or')
    end
end
ylim([0 1])
ylabel('Earlies')
pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-summary-' date '.pdf']);
exportfig_print(gcf, pn, 'FileFormat', 'pdf');

%earlies and threshold timecourse
figure; 
c= strvcat('k', 'r');
start = 1;
for imouse = 1:length(pv.mouse_mat)
    plot_mat = [];
    ipos = 1;
    if length(mouse(imouse).task(itask).pos(ipos).ind)>0
        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).ind)
            if ~isnan(mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh)
                pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).expt(iexp).n_early./mouse(imouse).task(itask).pos(ipos).expt(iexp).n_all;
                plot_mat = [plot_mat; iexp, pct_earlies_block1, mouse(imouse).task(itask).pos(ipos).expt(iexp).thresh];
            end
        end
    end
    if length(plot_mat)>0
        subplot(3,2,start)
        [ax, h1, h2] = plotyy(plot_mat(:,1), plot_mat(:,2), plot_mat(:,1), plot_mat(:,3));
        set(h1, 'Color', c(1), 'Marker', 'o')
        set(h2, 'Color', c(2), 'Marker', 'o')
        set(ax, 'YLimMode', 'auto', 'YTickMode', 'auto', 'YTickLabelMode', 'auto', 'xlim', [1 size(plot_mat,1)])
        set(ax(2), 'YColor', c(2))
        set(ax(1), 'YColor', c(1))
        title(num2str(pv.mouse_mat(imouse)))
        start = 1+start;
    end
end
pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-tcourse-task' num2str(itask) '-' date '.pdf']);
exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end