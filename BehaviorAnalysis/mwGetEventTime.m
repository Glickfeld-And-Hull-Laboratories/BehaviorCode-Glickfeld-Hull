function time_us = mwGetEventTime(events, codec, tag, ...
                                   occurrence, value, ignoreMissingStr)
%mwGetEventTime: given tag and possibly value, return event time
%
%  time_us = subGetEventTime(events, codec, tag, ...
%                            occurrence, value, ignoreMissingStr)
%
%
%  MH 100115: created
%$Id$

if nargin < 4 || isempty(occurrence), occurrence = 'all'; end
if nargin < 5 || isempty(value), value = []; end
if nargin < 6, ignoreMissingStr = []; end

% check ignoreMissingStr
if isempty(ignoreMissingStr)
  doIgnoreMissing=false;
elseif (isnumeric(ignoreMissingStr) || islogical(ignoreMissingStr)) && ignoreMissingStr == 0 
  doIgnoreMissing=false; 
elseif (isnumeric(ignoreMissingStr) || islogical(ignoreMissingStr)) && ignoreMissingStr == 1
  doIgnoreMissing=true; 
elseif ischar(ignoreMissingStr) && strcmpi(ignoreMissingStr, 'ignoremissing')
  doIgnoreMissing=true;
else
  error('Unknown ignoremissing value: %s', mat2str(ignoreMissingStr));
end

% check value
if ischar(value)
    if ~strcmpi(value, 'any')
        error('value must be numeric or ''any''');
    end
end


codes = [events.event_code];
tCode = codec_tag2code(codec, tag);
eventNs = find(codes == tCode);

if length(eventNs) < 1
  if ~doIgnoreMissing
    disp(sprintf('Code not found (in mwGetEventTime): %s', tag));
  end
  time_us = []; 
  return
else
  % at least one
  eVals = cat(2,events(eventNs).data);
  eTs = cat(2,events(eventNs).time_us);

  % do indexing based on value if desired
  if isempty(value) || (ischar(value) && strcmpi(value) == 'any')
      dataNs = 1:length(eVals);
      valStr = 'any';
  else
      dataNs = find(eVals == value);
      valStr = char(value);
      if length(dataNs) == 0
          disp(sprintf('Code %s with value %s not found', tag, valStr));
          time_us = [];  
          return
      end
  end

  if strcmpi(occurrence, 'all') || isempty(occurrence)
      time_us = eTs(dataNs);
      return
  elseif strcmpi(occurrence, 'last')
      time_us = eTs(dataNs(end));
      return
  elseif length(dataNs) < occurrence
      disp(sprintf('Asked for %d codes %s with value %d, but found %d with this value (%d total)', ...
                   occurrence, tag, value, length(dataNs), length(eventNs)));
      time_us = [];
      return
  else
      time_us = eTs(dataNs(occurrence));
      return
  end
end


