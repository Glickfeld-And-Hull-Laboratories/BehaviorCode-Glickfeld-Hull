function assert_args (expr, errmsg)
%ASSERT_ARGS (ps-utils): test condition.  Use for checking arguments!
%   ASSERT_ARGS (EXPR, ERRMSG)
%
%   Calls ASSERT with errmsg set to 'PSutils:argumentError'
% 
%  MH - http://github.com/histed/tools-mh

errId = 'PSutils:argumentError';

if nargin < 2, errmsg = 'An input argument is invalid'; end
assert(expr, errmsg, errId, 'cannotdisable');
