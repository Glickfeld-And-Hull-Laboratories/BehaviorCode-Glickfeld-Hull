%% create your own script to do this if you want to export elsewhere or keep
%% figures around

%% run one at a time
while true
    close all;
    
    [b figH fitS bootS] = doSingleFitFromIndex(...
        'DoSubjAndDate', 'next', ...
        'NBootstrapReps', 1000, ...
        'LoadDataDebug', true, ...
        'DoExport', true);
    drawnow;
    if isempty(b), break; end
     
    figListH = get(0, 'Children');
    close(setdiff(figListH, figH));
end

