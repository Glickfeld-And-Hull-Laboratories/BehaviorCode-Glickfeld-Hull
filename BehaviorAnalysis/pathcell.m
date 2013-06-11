function cellP = pathcell
%PATHCELL (ps-utils): return current path in a cellstr
%   cellP = PATHCELL
%
%  MH - http://github.com/histed/tools-mh

p = path;
p1 = strrep(p, ':', ''',''');
cellP = eval(['{''', p1, '''}']);


