function outstr = strip_chars_except (instr, chars_to_not_strip);
%strip_chars_except (ps-utils): remove all chars not specified from str
%
% USAGE
%  function outstr = strip_chars_except (instr, chars_to_not_strip);
%
%   MH - http://github.com/histed/tools-mh

outstr = [];
for i = 1:size(instr,1)
  boolsleft = ismember (instr(i,:), chars_to_not_strip);
  indicesleft = find(boolsleft);
  outline = instr (i,indicesleft);
  outstr = strvcat (outstr, outline);
end

