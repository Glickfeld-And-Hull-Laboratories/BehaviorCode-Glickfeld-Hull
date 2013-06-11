function join_pdf_files(inputNameCellstr, outputName, doDeleteInput);
%JOIN_PDF_FILES(inputNameCellstr, outputName, doDeleteInput)
%
%   inputNameCellstr can be a pattern
%   uses an adhoc method on each platform
%
%   This is just a stub to later be expanded - works on the mac w/ pdftk
%   installed; better pattern/cellstr handling should be done
%
%  MH - http://github.com/histed/tools-mh

if nargin < 3, doDeleteInput = false; end

if ~iscell(inputNameCellstr)
    inputNamePat = inputNameCellstr; % assume it's a pattern
else
    inputNamePat = sprintf('%s ', inputNameCellstr{:});
end


if isunix
    [s,w] = unix(sprintf('/usr/local/bin/pdftk %s cat output %s 2>&1', ...
                 inputNamePat, outputName));
    if s ~= 0
        error('%s: pdftk error: result %d message %s', mfilename, s, w);
    end

    % above will raise error on problems, so this is safe
    if doDeleteInput
        [s,w] = unix(sprintf('rm %s', inputNamePat));
        if s2abcd~=0
            error('%s: rm error: result %d message %s', mfilename, s, w);
        end
    end
end

