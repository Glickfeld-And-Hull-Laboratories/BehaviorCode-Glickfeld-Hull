function plotEarlyRegress(itask)

    mouse = createExptStructPV;
    pv = behavParamsPV;
    rc = behavConstsHADC8;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);


    for imouse = 1:length(pv.mouse_mat)
        for ipos = 1:length(pv.pos_mat)
            for ipow = 1:length(pv.power_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
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
                    RTs = cell2mat(input.reactTimesMs);
                    block1_ind_early = find(RTs(block1_ind_all)<=125);
                    block2_ind_early = find(RTs(block2_ind_all)<=125);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_all = length(block1_ind_all);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_all = length(block2_ind_all);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).n_early = length(block1_ind_early);
                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).n_early = length(block2_ind_early);                    
                end
            end
        end
    end

    % percent earlies by threshold by mouse
    c = strvcat('k', 'r');
    figure;
    for iblock = 1:2
        start = 1;
        b_mat = zeros(3, length(pv.mouse_mat));
        p_mat = zeros(3,2,length(pv.mouse_mat));
        for imouse = 1:length(pv.mouse_mat)
            z =0;
            x_mat = [];
            y_mat = [];
            pow_mat = [];
            ipos = 1;
            for ipow = 1:length(pv.power_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)
                        pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                        z = 1;
                        y_mat = [y_mat; mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh];
                        x_mat = [x_mat; pct_earlies_block1 pv.power_mat(ipow) 1];
                    end
                end
            end
            if z == 1;
                logy_mat = log10(y_mat);
                subplot(3,2,start)
                plot(x_mat(:,1), logy_mat, ['o' c(iblock)])
                hold on
                [b, p] = regress(logy_mat,x_mat);
                r = b(1) .*x_mat(:,1) + b(3);
                plot(x_mat(:,1), r, ['-' c(iblock)]);
                title([num2str(pv.mouse_mat(imouse))])
                start = 1+start;
                ylabel('Thresh')
                xlabel('Percent Earlies')
                ylim([-1.2 0])
                xlim([0 1])
            end 
            if iblock == 2
                b_mat(:,imouse) = b;
                p_mat(:,:,imouse) = p;
            end
        end
    end
    pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-task' num2str(itask) '-regressions-allpowers-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');

    % percent earlies by threshold by mouse at change in thresh = ~2
    c = strvcat('k', 'r');
    figure;
    b_mat = zeros(2,length(pv.mouse_mat));
  	p_mat = zeros(2,2,length(pv.mouse_mat));
    for iblock = 1:2
        start = 1;
        for imouse = 1:length(pv.mouse_mat)
            if pv.chr2_mat(imouse) == 1
                z =0;
                x_mat = [];
                y_mat = [];
                pow_mat = [];
                ipos = 1;
                pows = [0.25 0 0.25 0.15 0];
                pow = pows(imouse);
                for ipow = find(pv.power_mat == pow)
                    for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                        if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)
                            pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                            z = 1;
                            y_mat = [y_mat; mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh];
                            x_mat = [x_mat; pct_earlies_block1];
                        end
                    end
                end
                if z == 1;
                    logy_mat = log10(y_mat);
                    subplot(3,2,start)
                    plot(x_mat, logy_mat, ['o' c(iblock)])
                    hold on
                    x_mat_1s = cat(2, x_mat, ones(length(x_mat),1));
                    [b, p] = regress(logy_mat, x_mat_1s);
                    r = b(1) .*x_mat + b(2);
                    plot(x_mat, r, ['-' c(iblock)]);
                    title([num2str(pv.mouse_mat(imouse)) ' power = ' num2str(pow) ' mW'])
                    text(0.1, (0.1-(iblock*.2)), ['R = ' num2str(chop(b(1),2)) 'x +' num2str(chop(b(2),2))], 'Color', c(iblock)); 
                    start = 1+start;
                    ylabel('Thresh')
                    xlabel('Percent Earlies')
                    xlim([0 1])
                    ylim([-1.2 0])
                    b_mat(iblock,imouse) = b(1);
                    p_mat(iblock,:,imouse) = p(1,:);
                end
            end
        end
    end
    pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-task' num2str(itask) '-regressions-onepower-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
    
    % slope of delta(earlies) vs delta(threshold) by mouse at change in thresh = ~2
    slope = [];
    for imouse = 1:length(pv.mouse_mat)
        if pv.chr2_mat(imouse) == 1
            slope_mat = [];
            x_mat = [];
            y_mat = [];
            ipos = 1;
            pows = [0.25 0 0.25 0.15 0];
            pow = pows(imouse);
            for ipow = find(pv.power_mat == pow)                    
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh)
                        if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh)
                            for iblock = 1:2
                                pct_earlies(iblock) = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                                thresh(iblock) = log10(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh);
                                if iblock == 1
                                    x_mat = [x_mat; pct_earlies(iblock) 1];
                                    y_mat = [y_mat; thresh(iblock)];
                                end
                            end
                            slope_mat = [slope_mat; (thresh(2)-thresh(1))./(pct_earlies(2)-pct_earlies(1))];
                        end
                    end
                end
            end
            slope(imouse).indiv = slope_mat;
            slope(imouse).indiv_avg = mean(slope_mat,1);
            slope(imouse).indiv_med = median(slope_mat,1);
            [b, p] = regress(y_mat,x_mat);
            slope(imouse).one_slope = b(1);
            slope(imouse).one_ci95 = p(1,:);
        end
    end
    
    
    % percent earlies by power 
    c = strvcat('k', 'r');
    figure;
    iblock = 1;
    start = 1;
    for imouse = 1:length(pv.mouse_mat)
        if pv.chr2_mat(imouse) == 1
            z =0;
            x_mat = [];
            y_mat = [];
            ipos = 1;
            for ipow = 1:length(pv.power_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)
                        pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                        z = 1;
                        x_mat = [x_mat; pct_earlies_block1 1];
                        y_mat = [y_mat; log10(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)];
                    end
                end
            end
            if z == 1;
                subplot(3,3,start)
                plot(x_mat(:,1), y_mat, ['o' c(iblock)])
                hold on
                [b, p] = regress(y_mat,x_mat);
                r = b(1) .*x_mat(:,1) + b(2);
                plot(x_mat(:,1), r, ['-' c(iblock)]);
                title([num2str(pv.mouse_mat(imouse))])
                text(.1, -.2, ['R= ' num2str(chop(b(1),2))]);
                ylabel('log(Thresh)')
                xlabel('FAs')
                ylim([-1.2 0])
                xlim([0 0.8])
                start = start +3;
                slope(imouse).all_slope = b(1);
                slope(imouse).all_ci95 = p(1,:);
            end
        end
    end
    
    iblock = 2;
    start = 2;
    for imouse = 1:length(pv.mouse_mat)
        if pv.chr2_mat(imouse) == 1
            z =0;
            x_mat = [];
            FA_mat = [];
            thresh_mat = [];
            ipos = 1;
            for ipow = 1:length(pv.power_mat)
                for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                    if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)
                        pct_earlies_block1 = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                        z = 1;
                        FA_mat = [FA_mat; pct_earlies_block1];
                        thresh_mat = [thresh_mat mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh];
                        x_mat = [x_mat; pv.power_mat(ipow) 1];
                    end
                end
            end
            if z == 1;
                subplot(3,3,start)
                y_mat = log10(thresh_mat)';
                plot(x_mat(:,1), y_mat, ['o' c(iblock)])
                hold on
                [b, p] = regress(y_mat,x_mat);
                r = b(1) .*x_mat(:,1) + b(2);
                plot(x_mat(:,1), r, ['-' c(iblock)]);
                text(.1, -.2, ['R= ' num2str(chop(b(1),2))]);
                ylabel('log(Thresh)')
                xlabel('Power')
                ylim([-1.2 0])
                xlim([0 .6])
                start = start +1;
                subplot(3,3,start)
                y_mat = FA_mat;
                plot(x_mat(:,1), y_mat, ['o' c(iblock)])
                hold on
                [b, p] = regress(y_mat,x_mat);
                r = b(1) .*x_mat(:,1) + b(2);
                plot(x_mat(:,1), r, ['-' c(iblock)]);
                text(.1, .8, ['R= ' num2str(chop(b(1),2))]);
                ylabel('FAs')
                xlabel('Power')
                ylim([0 1])
                xlim([0 .6])            
                start = 2+start;
            end
        end
    end
    pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-task' num2str(itask) '-regressions-power-thresh-FAs-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
    
    s = [];
    s.mouse_mat = [];
    s.pow_mat = [];
    s.FA_mat = [];
    s.thresh_mat = [];
    for imouse = 1:length(pv.mouse_mat)
        for ipow = 1:length(pv.power_mat)
            for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh)
                    s.mouse_mat = [s.mouse_mat; pv.mouse_mat(imouse)];
                    s.pow_mat = [s.pow_mat; pv.power_mat(ipow)];
                    FA_mat = zeros(1,2);
                    thresh_mat = zeros(1,2);
                    for iblock = 1:2
                        FA_mat(1,iblock) = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_early./mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).n_all;
                        thresh_mat(1,iblock) = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(iblock).thresh;
                    end
                    s.FA_mat = [s.FA_mat; FA_mat];
                    s.thresh_mat = [s.thresh_mat; thresh_mat];
                end
            end
        end
    end
    pn = fullfile(rc.fitOutputSummary, 'FA_thresh_power_mouse.mat');
    save(pn, 's')
    
    figure;
    for imouse = 1:length(pv.mouse_mat)
        if pv.chr2_mat(imouse) == 1
            subplot(2,2,1)
            plot(imouse, slope(imouse).one_slope,'ob')
            hold on
            errorbar(imouse, slope(imouse).one_slope, slope(imouse).one_ci95(1)-slope(imouse).one_slope, slope(imouse).one_slope-slope(imouse).one_ci95(2))
            hold on
            plot(imouse,slope(imouse).indiv_med, 'ok')
            ylim([-5 1])
            ylabel('slope')
            title('One power')
            subplot(2,2,2)
            plot(imouse, slope(imouse).all_slope,'ob')
            hold on
            errorbar(imouse, slope(imouse).all_slope, slope(imouse).all_ci95(1)-slope(imouse).all_slope, slope(imouse).all_slope-slope(imouse).all_ci95(2))
            hold on
            plot(imouse,slope(imouse).indiv_med, 'ok')
            ylim([-5 1])
            title('All powers')
        end
    end
    pn = fullfile(rc.fitOutputSummary, ['Earlies-thresh-task' num2str(itask) '-slope-diff-' date '.pdf']);
    exportfig_print(gcf, pn, 'FileFormat', 'pdf');
end