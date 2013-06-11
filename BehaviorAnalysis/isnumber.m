function outbool = isnumber (instr)
%isnumber (ps-utils): true for numbers in strings.
%
% USAGE
%  function outbool = isnumber ()
%
% true for characters like '1', false for [1].
% See also ISLETTER.
% This is likely slow, isletter is implemented as an internal
% function.
%
%   MH - http://github.com/histed/tools-mh

numbers = '1234567890';
outbool = [];

if ~ischar(instr)
  warning ('Parameter is not a string, returning []');
  return
end

outbool = ismember (instr, numbers);



