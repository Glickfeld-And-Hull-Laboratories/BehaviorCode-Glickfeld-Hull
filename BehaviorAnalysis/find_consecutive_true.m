function [starti,endi] = find_consecutive_true (invec, how_many)
%FIND_CONSECUTIVE_TRUE (ps-utils): indices of consecutive true entries
%  [STARTI ENDI] = FIND_CONSECUTIVE (INVEC, HOW_MANY)
%
%  find consecutive true entries of HOW_MANY or more in a logical vector
%
%    MH - http://github.com/histed/tools-mh

if nargin < 2, how_many = 2; end

assert(islogical(invec));

[startS endS] = find_consecutive(invec, 0);

runLens = endS-startS;
isTrueToStart = invec(startS) == true;

desIx = runLens >= how_many & isTrueToStart;

starti = startS(desIx);
endi = endS(desIx);





