function [] = plotReqXHolds( ds )
% Used to plot x = hold times y = required hold time when input is a
% behavior structure for HADC8.xml
disp(isfield(ds, 'randReqHoldMaxMs'))
if isfield(ds, 'randReqHoldMaxMs'),
    hts = cell2mat_padded(ds.holdTimesMs);
    req = cell2mat_padded(ds.tTotalReqHoldTimeMs)+ds.tooFastTimeMs;
    
    [rewOnset, Ix] = sort(req);
    rewEnd = cell2mat_padded(ds.tTotalReqHoldTimeMs)+ds.reactTimeMs;
    rewEnd = rewEnd(Ix);    
    
    [a] = scatter(req, hts);
    jbfill(rewOnset, rewEnd, rewOnset, 'g', 'g', 1, 0.25)
    ylabel('Hold Time (ms)');
    xlabel('Required Hold Time (ms)');
    grid on
    title(ds.savedDataName(40:end));
    outPath = '~/Desktop/HoldXReq/';
    fName = strcat('reqXholds-', ds.savedDataName(40:55), '.pdf')
    fileName = strcat(outPath, fName);
    exportfig_print(gcf, fileName, 'FileFormat', 'pdf', ...
             'Size', [12 9], ...
             'PrintUI', false)
    pdf = gcf;
end
end