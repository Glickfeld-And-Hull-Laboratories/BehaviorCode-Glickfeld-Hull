function frm_javasetup(poiPath);
%FRM_JAVASETUP: configure java path for xls manipulation
%
%   frm_javasetup(poiPath)
%  
%   example:
%      frm_javasetup('~/path/to/poi-3.8')
%
%   This will automatically add the required jar files to the java path, including the Apache POI library 
%     (used for lowlevel XLS/XLSX reading and writing) and frm_ java subroutines called by the .m files.
%
%  MH - http://github.com/histed/tools-mh

% deal with bug parsing '~' dirs on Mac
if any(poiPath == '~')
    [~,b0] = fileattrib(poiPath);
    poiPath = b0.Name;
end

ds = dir(fullfile(poiPath, 'poi-*.jar'));
names = {ds.name};
if length(names) < 3
    error('POI jar files not found; you gave path %s', poiPath);
end


jarPats = { ...
    'poi-[0-9\.-]*.jar', ...
    'poi-ooxml-[0-9\.-]*.jar', ...
    'poi-ooxml-schemas-[0-9\.-]*.jar' };

for iR = 1:3  
    % 3 total jar files to find here
    tPat = jarPats{iR};
    nameR = regexpi(names, tPat);
    matchIx = ~cellfun(@isempty, nameR);
    assert(sum(matchIx) == 1, 'bug: jar file not found');
    subJavaAddIfMissing(fullfile(poiPath, names{matchIx}));
end

% manually add one in another dir
d2 = fullfile(poiPath, 'ooxml-lib');
n0 = dir(fullfile(d2, 'xmlbeans-*.jar'));
xbName = {n0.name};
assert(length(xbName) == 1, 'bug: xmlbeans jar not found');
subJavaAddIfMissing(fullfile(d2, xbName{1}));

% add the java stub helper function
fullNameToThis = which(mfilename);
tPath = fileparts(fullNameToThis);
subJavaAddIfMissing(fullfile(tPath, 'java/frm_jcReadXls'));

function subJavaAddIfMissing(tPath)
jp = javaclasspath;
if ~any(strcmp(jp, tPath))
    javaaddpath(tPath);
end
