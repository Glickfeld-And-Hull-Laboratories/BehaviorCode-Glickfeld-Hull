function boolOut = isstropt(putative_stropt)
%ISSTROPT (ps-utils): does this parameter look like a stropt cell array?
%   boolOut = ISSTROPT(PUTATIVE_STROPT)
%   Used for stropt param checking
%
%  MH - http://github.com/histed/tools-mh

if ~(isempty(putative_stropt) ...
     || (mod(length(putative_stropt),2) == 0 ...
         && all(cellfun('isclass',putative_stropt(1:2:end),'char'))))
    boolOut = false;
else
    boolOut = true;
end
