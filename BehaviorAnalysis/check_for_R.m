function check_for_R 
%check_for_R (ps-utils): See if R is initialized, if not, set it up
%
%   Uses the RMatlab interface.  If you call initializeR twice in the same
%   Matlab session, matlab crashes.
%
%  MH - http://github.com/histed/tools-mh

global RMATLAB_R_IS_INITIALIZED

if isempty(RMATLAB_R_IS_INITIALIZED)
    % not started yet, start it
    
    retval = initializeR({'Rmatlab' '--no-save', '--silent', '--no-restore'});
    if retval ~= 1
        error('Initializing R: retval %d', retval);  
    end
    
    RMATLAB_R_IS_INITIALIZED = true;
else
    % already started, just return
end

    
    
