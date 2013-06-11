function assert_always(varargin)
%ASSERT_ALWAYS (ps-utils): test condition, if it fails, error (cannot disable)
%  ASSERT_ALWAYS(EXPR, ERRMSG, ERRID) 
%
%  Setting the global NO_ASSERTIONS variable has no effect on this function
%
%  See also ASSERT
% 
%  MH - http://github.com/histed/tools-mh

if nargin < 3, varargin{3} = 'PSUTILS:always_assert_error'; end
% this handles the case where nargin == 1 because setting the 3rd element
% makes the 2d empty if it doesn't already exist

assert(varargin{1:3}, 'cannotdisable');