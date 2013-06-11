function merged = struct_union(primary, secondary)
%STRUCT_UNION (ps-utils): merge two structures into one.
%   MERGED = STRUCT_UNION(PRIMARY, SECONDARY)
%   primary, secondary: structures to be merged.  If a field exists in both
%   structs, the value from primary will be used.
%
%   Notes: could be optimized for speed, now O(n) in the number of fields in
%   primary.
%
%  MH - http://github.com/histed/tools-mh

if isempty(secondary), merged=primary; return; end
if isempty(primary), merged=secondary; return; end
merged = secondary;
s = fieldnames(primary);
for i = 1:length(s)
   merged = setfield (merged, s{i}, getfield(primary, s{i}));
end
