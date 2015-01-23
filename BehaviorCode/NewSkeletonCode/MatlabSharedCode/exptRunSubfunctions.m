function input = exptRunSubfunctions(data_struct, input, subFHCell)
% call at end of main experiment Matlab script
%
% each function handle in subFHCell is called as:
%     input = fh(ds, input)
%
% MH 130107

%% consts
doProfile = false;

%%

if doProfile
    profile on;
end

tic
try
    if ~(input.experimentXmlTrialId==15), %completely skip 
        input = exptSaveMatlabState(data_struct, input);
    else
        disp('Lego Detected, must save manually')
    end

    % do the plot
    if ~isfield(input, 'lastPlotElapsedS')
        input.lastPlotElapsedS = 0;
    end
    if input.lastPlotElapsedS > 2
        % skip every other plot to make sure we catch up
        disp(sprintf('Skipping plot, last elapsed was %5.1fs', ...
                     input.lastPlotElapsedS));
        input.lastPlotElapsedS = 0;
    else
        % last one was fast, plot this time
        tic
        nFH = length(subFHCell);
        for iF = 1:nFH
            % this input is the CHANGED event stream w/ first trial consts extracted
            input = feval(subFHCell{iF}, data_struct, input);
        end
        
        input.lastPlotElapsedS = toc;
        %disp(input.lastPlotElapsedS);
    end
    
catch ex
    disp('??? Error in subfunction; still saving variables for next trial')
    printErrorStack(ex);
end
tEl = toc;
%if tEl > 1.3
%  disp(sprintf('          Elapsed: %gsec  ', chop(tEl,2)));
%end

if doProfile
    p = profile('info');
    save('/tmp/profileHADC8.mat', p)
end
