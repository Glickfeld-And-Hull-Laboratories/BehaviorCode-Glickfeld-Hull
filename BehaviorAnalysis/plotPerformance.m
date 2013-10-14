function [x] = plotPerformance(subjs)
global dataStru
subPlotSz = {2,1};
for i=1:length(subjs),
    sub = subjs(i);
    pCorr = dataStru(sub).values.perCorrects;
    pEarly = dataStru(sub).values.perEarlies;
    pLate = dataStru(sub).values.perLates;
    pNR = dataStru(sub).values.perNR;   
    medH = dataStru(sub).values.medianHold;
    winS = dataStru(sub).values.windowStart;
    uSTD = dataStru(sub).values.uSTD;
    lSTD = dataStru(sub).values.lSTD;
    nDays = length(medH);
    xDays = 1:nDays;
    
    axH = subplot(subPlotSz{:}, 1);
    
    hold on
    plot(pCorr, 'k', 'LineWidth', 2);
    plot(pEarly, 'c', 'LineWidth', 2);
    plot(pLate, 'm', 'LineWidth', 2);
    plot(pNR, 'Color', [1 0.64 0], 'LineWidth', 2);
    ylabel('Percent');
    xlabel('Training Day');
    tName = strcat('Trial Outcome Performance -- i', num2str(sub), '-- Generated:', datestr(today, 'dd mmmm yyyy'));
    title(tName);
    axis tight
    ylim([0 1]);
    legend('Correct', 'Early', 'Late', 'No Release') 
    grid on

    axH = subplot(subPlotSz{:}, 2);
    
    SD = jbfill(xDays, uSTD, lSTD, [0 0 1], [0 0 1], 1, 0.5);
    %plot(uSTD, 'b--', 'LineWidth', 1)
    %plot(lSTD, 'b--', 'LineWidth', 1)
    hold on
    medPlot = plot(medH, 'k', 'LineWidth', 2);
    winPlot = plot(winS, 'r', 'LineWidth', 2);
    ylabel('Hold Time (ms)');
    xlabel('Training Day');
    axis tight
    grid on
    ylim([0 2500]);
    tName = strcat('Lever Hold Performance -- i', num2str(sub));%' -- Generated:', datestr(today, 'dd mmmm yyyy'));
    title(tName);
    legend('\pm 1 SD', 'Median Hold Time', 'Reward Window Start');

    sName = strcat('graphPerformance-i', num2str(sub),'-', datestr(today, 'yymmdd'), '.pdf');
    saveas(gcf, sName)
    close
end