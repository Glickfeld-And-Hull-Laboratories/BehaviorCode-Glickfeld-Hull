function dsOut = convertDataSerial(ns)
% Sanitize the serial data structure to be the same as sanitized matlab
% pass in the 'ns' struct
% MH 130118

td = ns.parsedTrialData;
ds = td.endStateS;


% backward compat
ds.holdTimesMs = ds.actualHoldTimeMs;
ds.reactTimesMs = ds.reactTimeMs;

% process dates
ds.startDateVec = ns.MetaTags.DateTimeRaw(1:5);
dr = ns.MetaTags.DateTimeRaw;
ds.startDateVec(6) = dr(6) + dr(7)/1e3 + dr(8)/1e6;

% compute outcomecell
ds.trialOutcomeCell = subComputeTrialOutcomeCell(ds);

dsOut = ds;
