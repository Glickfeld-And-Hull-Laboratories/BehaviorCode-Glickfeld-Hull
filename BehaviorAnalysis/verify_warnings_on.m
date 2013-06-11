function verify_warnings_on
%VERIFY_WARNINGS_ON (ps-utils): if warning state is off, error()
%   VERIFY_WARNINGS_ON
%
%  MH - http://github.com/histed/tools-mh

ws = warning;
wId = { ws.identifier };

allState = ws(strcmp(wId, 'all')).state;

if strcmp(allState, 'off')
    error('Global warning state is off');
else
    assert(strcmp(allState, 'on') || strcmp(allState, 'backtrace'), ...
           'Warning state is not ''off'' ''on'' or ''backtrace''');
end


