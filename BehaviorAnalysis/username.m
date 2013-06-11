function outName = username
%USERNAME (ps-utils): return user name of the caller
%   T = USERNAME returns user name of the caller 
%
% by Vincent Bonin 
% based on HOSTNAME by Mark Histed
%  MH - http://github.com/histed/tools-mh

outName = getenv('USERNAME');

return;

% PREVIOUS CODE, DID NOT WORK ON XP PRO HOSTS.

% % cache check
% [retval,tName] = system('whoami');
% 
% % some error checks
% assert(retval == 0, 'Error running hostname');
% 
% [pathstr,name]=fileparts(tName);
% 
% outName = deblank(name);
