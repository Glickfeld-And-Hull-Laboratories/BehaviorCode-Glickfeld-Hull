function dirTime = lastmodtime(filename)
%LASTMODTIME (ps-utils): Use DIR to get the last modification time of a file
%
%   dirTime = LASTMODTIME(filename)
%
%  MH - http://github.com/histed/tools-mh

d = dir(filename);
if isempty(d)
    error('File not found');
else
    dirTime = d.date;
end

