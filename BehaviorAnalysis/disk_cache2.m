function outD = disk_cache2(opStr, varargin)
% DISK_CACHE2 (mh-tools) cache data to disk as identifiable mat files, specify directory
%
%   diskCache2('get', cacheKeyStr, cacheDir, allowMissing)
%   diskCache2('set', data, cacheKeyStr, cacheDir, quietStr)
%
%   You are responsible for creating cacheKeyStr and making sure it's unique
%
% histed 111005


% arg processing based on opstr
switch lower(opStr)
  case 'set'
    argC = varargin;
    nDesArgs = 4;
    if length(varargin) < nDesArgs, varargin{nDesArgs} = []; end % pad
    [inData, cacheKeyStr, cacheDir, quietStr] = deal(varargin{1:nDesArgs});

    % complete empty
    if (length(quietStr)==1 && quietStr == true) || strcmpi(quietStr, 'quiet')
        quiet = true;
    elseif (length(quietStr)==1 && quietStr == false) || isempty(quietStr)
        quiet = false;
    else 
        error('bad value for quietStr: %s', mat2str(quietStr));
    end

  
  case 'get'
    argC = varargin;
    nDesArgs = 3;
    if length(varargin) < nDesArgs, varargin{nDesArgs} = []; end % pad
    [cacheKeyStr, cacheDir, allowMissingStr] = deal(varargin{1:nDesArgs});

    % complete any empty arguments
    if (length(allowMissingStr)==1 && allowMissingStr == true) || strcmpi(allowMissingStr, 'allowMissing')
        allowMissing = true;
    elseif (length(allowMissingStr)==1 && allowMissingStr == false) || isempty(allowMissingStr)
        allowMissing = false;
    else 
        error('bad value for allowMissingStr: %s', mat2str(allowMissingStr));
    end
  
  otherwise
    error(sprintf('invalid opStr %s', opStr));
end

if isempty(cacheDir)
    cacheDir = getTmpfile;
end



% use keystr directly
outFileName = fullfile(cacheDir, ['diskcache2_' cacheKeyStr '.mat']);

switch lower(opStr)
  case 'set'
    save(outFileName, 'inData', '-v7.3');  % v7.3, support >2GB
    dirS = dir(outFileName);
    w = whos('inData');
    if ~quiet
        fprintf(1, 'disk_cache2: saved %sb (%sb compressed), file %s\n', ...
                num2str_metric(w.bytes,2), ...
                num2str_metric(dirS.bytes,2), ...
                outFileName);
    end
    
  case 'get'
    if ~exist(outFileName, 'file')
        if allowMissing
            outD = [];
            return
        else
            error('Cache entry not found for key %s', cacheKeyStr);
        end
    end
    
    b = load(outFileName);
    outD = b.inData;
end

