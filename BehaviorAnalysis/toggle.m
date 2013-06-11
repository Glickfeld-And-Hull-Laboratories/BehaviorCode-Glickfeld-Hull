function new = toggle (old)
%TOGGLE (ps-utils): changes the state of a boolean var.
%   NEW = TOGGLE (OLD)
%   Examples:
%   toggle (1) == 0 
%   toggle (0) == 1
%   toggle (548) == 0
%
%  MH - http://github.com/histed/tools-mh

if old == 0
  new = 1;
else
  new = 0;
end

  