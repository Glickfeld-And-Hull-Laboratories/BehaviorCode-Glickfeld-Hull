function [x] = plotPerformance(subjs)
global dataStru
subPlotSz = {3,1};
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
    nTrials = dataStru(sub).values.nTrials;
    nCorr = dataStru(sub).values.nCorrects;
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
    legend('Correct', 'Early', 'Late', 'No Release') ;
    grid on;

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
    
    axH = subplot(subPlotSz{:}, 3);
    [Ax, H1, H2] = plotyy(xDays, nCorr, xDays, nTrials);
    set(H1,'LineWidth', 2);
    set(H1,'Color', 'k');
    set(H2,'LineWidth', 2);
    set(H2,'Color', 'g');
    set(Ax(1), 'XLim', [1 max(xDays)]);
    set(Ax(2), 'XLim', [1 max(xDays)]);
    set(Ax(2), 'YLim', [0 max(nTrials)]);
    set(Ax(1), 'YLim', [0 max(nTrials)]);
    set(Ax(1),'YColor', [0 0 0])
    xlabel('Training Day');
    set(get(Ax(1),'YLabel'), 'String', 'Correct Trials');
    set(get(Ax(2),'YLabel'), 'String', 'Total Trials');
    set(Ax(2), 'YTick', [0:200:ceil(max(nTrials)/200)*200]);
    set(Ax(1), 'YTick', [0:200:ceil(max(nTrials)/200)*200]);
    grid on
    title('Number of Correct and Total Trials')
    fName = strcat('Performance -- i', num2str(sub), '-- Generated:', datestr(today, 'dd mmmm yyyy'));
    set(gcf, 'Name', fName);
    
    
    sName = strcat('graphPerformance-i', num2str(sub),'-', datestr(today, 'yymmdd'), '.pdf');
    saveas(gcf, sName)
    close
end