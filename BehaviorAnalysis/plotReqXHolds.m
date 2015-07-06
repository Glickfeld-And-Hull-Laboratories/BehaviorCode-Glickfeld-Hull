%function [ds] = plotReqXHolds( ds )
% Used to plot x = hold times y = required hold time when input is a
% behavior structure for HADC8.xml
if ispc,
    disp('I AM A PC')
    addpath('C:\Users\jake\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis\');
elseif isunix,
    disp('I AM A MAC OR LINUX COMPUTER')
    addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
end
disp(isfield(ds, 'randReqHoldMaxMs'))
if isfield(ds, 'randReqHoldMaxMs'),
    hts = cell2mat_padded(ds.holdTimesMs);
    hts = double(hts);
    req = double(cell2mat_padded(ds.tTotalReqHoldTimeMs))+double(ds.tooFastTimeMs);
    
    [rewOnset, Ix] = sort(req);
    rewEnd = cell2mat_padded(ds.tTotalReqHoldTimeMs)+ds.reactTimeMs;
    rewEnd = rewEnd(Ix);
    corrIx = hts>req;
    
    subplotSize = {1,2};
    axH = subplot(subplotSize{:}, 1)
    eHT = hts(~corrIx);
    totalHTs = length(eHT);
    [a b] = hist(eHT, ceil(range(eHT)/10));
    probs = a/totalHTs;
    cdfStuff = cumsum(probs);
    plot(cdfStuff,b)
    
    %pH1 = cdfplot(hts(~corrIx));
    %set(gca,'YTickLabel',[]);
    %set(gca,'XTickLabel',[]);
    xlabel('F(x)')
    ylabel('Hold Time of Early Releases(ms)')
    axis tight
    grid on
    ylim([min(hts) max(hts)])
    axH = subplot(subplotSize{:}, 2)
    pH2 = scatter(req(~corrIx), hts(~corrIx), 'b');
    %[a] = scatterhist(req(~corrIx), hts(~corrIx), 'Direction', 'out', 'NBins', [20 ceil(range(hts(~corrIx))/100)], 'Color', 'b');
    hold on
    pH3 = scatter(req(corrIx), hts(corrIx), 'k');
    if range(req)>50,
        xlim([min(req) max(req)]);
    end
    ylim([min(hts) max(hts)]);
    jbfill(rewOnset, rewEnd, rewOnset, 'g', 'g', 1, 0.25)
    ylabel('Hold Time (ms)');
    xlabel('Required Hold Time (ms)');
    grid on
    title(ds.savedDataName(40:end));
    outPath = 'C:\Users\jake\TempData\';
    fName = strcat('reqXholds2-', ds.savedDataName(40:55), '.pdf')
    fileName = strcat(outPath, fName);
    exportfig_print(gcf, fileName, 'FileFormat', 'pdf', ...
             'Size', [12 9], ...
             'PrintUI', false)
    pdf = gcf;
end
%end