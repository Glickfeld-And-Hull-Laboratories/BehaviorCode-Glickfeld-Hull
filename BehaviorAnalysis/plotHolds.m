sub = [001]
    
    %meanH = dataStru(sub).values.meanHold;
    medH = dataStru(sub).values.medianHold;
    uSTD = dataStru(sub).values.uSTD;
    lSTD = dataStru(sub).values.lSTD;
nDays = length(medH);
xDays = 1:nDays;
    
    hold on
    plot(medH, 'k', 'LineWidth', 2)
    
    jbfill(xDays, uSTD, lSTD, [0 0 1], [0 0 1], 1, 0.5)
    %plot(uSTD, 'b--', 'LineWidth', 1)
    %plot(lSTD, 'b--', 'LineWidth', 1)
    
    ylabel('Hold Time (ms)')
    xlabel('Training Day')
    axis tight
    grid on