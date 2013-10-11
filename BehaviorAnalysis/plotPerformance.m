
sub = [001]
    
    pCorr = dataStru(sub).values.perCorrects
    pEarly = dataStru(sub).values.perEarlies
    pLate = dataStru(sub).values.perLates
    pNR = dataStru(sub).values.perNR
    
    hold on
    plot(pCorr, 'k', 'LineWidth', 2)
    plot(pEarly, 'c', 'LineWidth', 2)
    plot(pLate, 'm', 'LineWidth', 2)
    plot(pNR, 'Color', [1 0.64 0], 'LineWidth', 2)
    ylabel('Percent')
    xlabel('Training Day')
    axis tight
    ylim([0 1])
    grid on
    