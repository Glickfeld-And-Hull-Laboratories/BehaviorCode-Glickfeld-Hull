function wc = wildcard
%WILDCARD (ps-utils): return wildcard character for present platform
%
%  MH - http://github.com/histed/tools-mh

if isunix
  wc = '*';
else % assume windows!
  wc = '*.*';
end
  
 