function padmat = vec2padded(invec, startix, whichcols, endix, pad, ...
                             allowNanEntries)
%VEC2PADDED (ps-utils): convert ragged mat from vector to mat
%   PADMAT = VEC2PADDED(INVEC, STARTIX, WHICHCOLS, ENDIX, PAD, ALLOWNANENTRIES)
%   This function takes a ragged matrix stored as a vector and a set of start
%   indices (STARTIX) and converts it to a rectangular matrix, padding as needed.
%   WHICHCOLS specifies a subset of columns to extract and return (0 or empty
%       means use all)
%   ENDIX is the set of end indices, if not specified or empty, it's computed
%     from STARTIX
%     STARTIX and ENDIX can overlap.
%   PAD is the double to pad with, default is NaN.
%   ALLOWNANENTRIES if true means if startix is NaN, the corresponding row in
%   the output will be all NaNs.  If false (default) NaN input gives error
%
%   PADMAT is in  column-major form.
%
%   See also: VEC2PADBYCOL
%
%  MH - http://github.com/histed/tools-mh

if length(startix) == 1 
    padmat = colvect(invec(startix:endix));
    return
end
if nargin < 6, allowNanEntries = false; end
if nargin < 5, pad = NaN; end
if nargin < 4 | isempty(endix)
    % compute endix, use last element in invec as end of last vector
    endix = startix(2:end)-1;
    endix = [ endix(:); length(invec) ];
end
if (nargin < 3 || isempty(whichcols) ...
    || (length(whichcols) == 1 && whichcols == 0))
    whichcols = 1:length(startix); 
end
assert(max(whichcols) <= length(startix), ...
       'Not enough trials in start/endix for some vals of whichcols');

% work with column vectors
endix = colvect(endix); 
startix = colvect(startix);
rowlens = endix-startix+1;  
% select elements from startix to endix, INCLUSIVE.  So if startix == 1 and
% endix == 10, rowlens = 10-1+1 = 10, and we take elements 1 through 10
% inclusive, a total of 10 elements!!!


nrows = max(rowlens(whichcols));
ncols = length(whichcols);
padmat = zeros(nrows, ncols) .* pad;

% could optimize this with sub2ind / cumsum, but JIT compiler handles it for
% us in 6.5 and above.  Profiling shows that this code not a problem. (5/9/03)
invec = colvect(invec);

doCols = 1:ncols;
nanIx = isnan(startix) | isnan(endix);
if any(nanIx)
    if ~allowNanEntries
        error('NaNs in input');
    else
        doCols = find(~nanIx);
    end
end

for i = doCols(:)'
    c = whichcols(i);
    padmat(1:rowlens(c),i) = invec(startix(c):endix(c));
end

