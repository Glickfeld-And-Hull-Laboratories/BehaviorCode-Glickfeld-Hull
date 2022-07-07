function selectivity_v_stimvars_for_summary_figs(smooth_select, tdays, stim_info)

% [smooth_select, smooth_MeanDT, tdays, smooth_select2, tdays2, stim_info] = getTrainingHistory2(behav_path, mouse);


    %% Subplot
    figure;
    
    % Level
    subplot(2,2,1)
    LVL = [stim_info.LVL1; stim_info.LVL2*2; stim_info.LVL3*3; stim_info.LVL4*4];
    norm_total_LVL = rescale(sum(LVL,1));
    
    diff_idx_LVL = find(diff(norm_total_LVL) ~= 0) + 1;
    plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    hold on
    ylim([-0.1 1])
    hline(0.9)
    hline(0)
    xlabel('Training Days')
    ylabel('Selectivity')
    title('Total Performance History (MA)')
    yyaxis right
    plot(tdays(diff_idx_LVL), norm_total_LVL(diff_idx_LVL), '*r', 'LineWidth', 1.24);
    if ~isempty(diff_idx_LVL)
        vline(diff_idx_LVL, ':r')
    end
    ylabel('Stim Difficulty')
    
    % Size
    
    subplot(2,2,2)
    size = stim_info.size;
    diff_idx_size = find(diff(size) ~= 0) + 1;
    
    plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    hold on
    ylim([-0.1 1])
    hline(0.9)
    hline(0)
    xlabel('Training Days')
    ylabel('Selectivity')
    title('Total Performance History (MA)')
    yyaxis right
    ylabel('Size')
    set(gca,'YDir', 'reverse')
    plot(tdays(diff_idx_size), size(diff_idx_size), 'ob', 'LineWidth', 1.24);
    ylim([30 50])
    if ~isempty(diff_idx_size)
        vline(diff_idx_size, ":b")
    end
    
    
    % ECC
    
    subplot(2,2,3)
    eccentricity = stim_info.eccentricity;
    diff_idx_ecc = find(diff(eccentricity) ~= 0) + 1;
    
    plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    hold on
    ylim([-0.1 1])
    hline(0.9)
    hline(0)
    xlabel('Training Days')
    ylabel('Selectivity')
    title('Total Performance History (MA)')
    yyaxis right
    plot(tdays(diff_idx_ecc), eccentricity(diff_idx_ecc), '>m', 'LineWidth', 1.24);
    ylim([0 30])
    if ~isempty(diff_idx_ecc)
        vline(diff_idx_ecc, ':m')
    end
    ylabel('Eccentricity')
    
    % FB
    
    subplot(2,2,4)
    FB = stim_info.FB;
    diff_idx_FB = find(diff(FB) ~= 0) + 1;
    
    plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    hold on
    ylim([-0.1 1])
    hline(0.9)
    hline(0)
    xlabel('Training Days')
    ylabel('Selectivity')
    title('Total Performance History (MA)')
    yyaxis right
    plot(tdays(diff_idx_FB), FB(diff_idx_FB), '+', 'LineWidth', 1.24);
    set(gca,'YDir', 'reverse')
    if ~isempty(diff_idx_FB)
        vline(diff_idx_FB)
    end
    ylabel('FB ON/OFF')
    
    
    %% Single plot
    % eccentricity = rescale(stim_info.eccentricity);
    % diff_idx_ecc = find(diff(eccentricity) ~= 0) + 1;
    % 
    % size = rescale(stim_info.size);
    % diff_idx_size = find(diff(size) ~= 0) + 1;
    % 
    % 
    % figure;
    % plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    % hold on
    % ylim([-0.1 1])
    % hline(0.9)
    % hline(0)
    % xlabel('Training Days')
    % ylabel('Selectivity')
    % title('Total Performance History (MA)')
    % yyaxis right
    % a = plot(tdays(diff_idx_LVL), norm_total_LVL(diff_idx_LVL), '*r', 'LineWidth', 1.24, 'DisplayName', 'Difficulty');
    % vline(diff_idx_LVL, ':r')
    % b = plot(tdays(diff_idx_size), size(diff_idx_size), 'ob', 'LineWidth', 1.24, 'DisplayName', 'Size');
    % vline(diff_idx_size, ':b')
    % c = plot(tdays(diff_idx_ecc), eccentricity(diff_idx_ecc), '>m', 'LineWidth', 1.24, 'DisplayName', 'Eccentricity');
    % vline(diff_idx_ecc, ':m')
    % ylabel('Stim Variables (Scaled)')
    % lgd = legend([a, b, c], 'Location', 'bestoutside');

end