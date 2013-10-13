function [x] = plotHolds(subjs)
global dataStru
for i=1:length(subjs),
    sub = subjs(i);
    %meanH = dataStru(sub).values.meanHold;
    medH = dataStru(sub).values.medianHold;
    winS = dataStru(sub).values.windowStart;
    uSTD = dataStru(sub).values.uSTD;
    lSTD = dataStru(sub).values.lSTD;
    nDays = length(medH);
    xDays = 1:nDays;
    
    hold on
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
    tName = strcat('Lever Hold Performance -- i', num2str(sub), ' -- Generated:', datestr(today, 'dd mmmm yyyy'));
    title(tName);
    legend('\pm 1 SD', 'Median Hold Time', 'Reward Window Start');
    sName = strcat('graphHolds-i', num2str(sub),'-', datestr(today, 'yymmdd'), '.pdf');
    saveas(gcf, sName)
    close
end