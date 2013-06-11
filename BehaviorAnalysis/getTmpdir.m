function tmpdir = getTmpdir
% getTmpdir (mh-tools): crossplatform way to find temp directory ($TMPDIR)
%
% histed 11105

if isunix
    tmpdir = getenv('TMPDIR')
    if isempty(tmpdir)
        tmpdir = '/tmp';
    end
else
    % windows
    tmpdir = getenv('TEMP')
    if isempty(tmpdir)
        error('Temp directory not found - fix mfile');
    end
end

