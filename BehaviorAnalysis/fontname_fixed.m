function fixed_name = fontname_fixed
%FONTNAME_FIXED (ps-utils): get the name of a fixed-width font on this system
%  NAME = FONTNAME_FIXED
%
%  MH - http://github.com/histed/tools-mh

l = listfonts;

% must be lowercase
possFontNames = { 'lucidatypewriter', ... 
                  'courier new', ...
                  'courier' };

lowerL = lower(l);
for i=1:length(possFontNames)
    tIx = strcmp(lowerL, possFontNames{i});
    if any(tIx)
        desN = find(tIx);
        assert(length(desN) == 1);
        break;
    end
end
fixed_name = l{desN};