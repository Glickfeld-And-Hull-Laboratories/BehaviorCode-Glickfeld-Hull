function [bs ds x1d fitX1d xd fitXd fitS bootS] = getAllBehavDataUsingIndex(varargin)
% This function reads behavioral data, finding and filtering based on index information.
%
%
%  bs is from getBehavDataForDay - basic transforms and aggregation of behav data
%  ds: raw behav data from matlab
%  x1d, fitX2d: relevant rows from the index and fit spreadsheet
%  xd, fitXd: full spreadsheets
%  fitS: fit parameters from weibullFitWithBootstrap
%  bootS: bootstrap results each rep from weibullFitWithBootstrap
%
% histed 121003

%% arg processing
userDefs = { ...
    'Subject', [], ...
    'DateStr', '', ...
    'UniqueId', '', ...
    'LoadDataDebug', false, ...
    'FitXd', []};
    
uo = stropt2struct(stropt_defaults(userDefs, varargin));

%% get xls data and find params
rc = behavConstsHADC8;

xd = readIndexFile;

if isempty(uo.FitXd)
    fitXd = readFitOutputFile;
else
    fitXd = uo.FitXd;
end

%% find index row
xdN = frm_findrownum(xd, ...
    { 'Subject', uo.Subject, ...
    'DateStr', uo.DateStr, ...
    'UniqueId', uo.UniqueId, ...
    });
x1d = frm_extractrow(xd, xdN);

% fit row
fitRowN = frm_findrownum(fitXd, ...
    { 'Subject', uo.Subject, ...
    'DateStr', uo.DateStr, ...
    'UniqueId', uo.UniqueId }, ...
    'IgnoreMissing', true);
if isnan(fitRowN) 
    fitX1d = [];
else
    fitX1d = frm_extractrow(fitXd, fitRowN);
end

if isnan(x1d.DataBlock), x1d.DataBlock = 'max'; end  % default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse whichtrials
whichTrs = x1d.TrialRangeToUse;
if isempty(whichTrs) || all(isnan(whichTrs))
    whichTrs = [];
elseif ischar(whichTrs)
    try
        whichTrs = eval(whichTrs);
    catch mexc
        warning('Error parsing text WhichTrials specification: ''%s'' - matlab error follows:', whichTrs);
        rethrow(mexc);
    end
end

%% parse other params
if isempty(x1d.MergeBlock1And2)
    mergeB1And2 = false;
else
    mergeB1And2 = x1d.MergeBlock1And2;
end

%% transforms
xT1FH = [];
xT2FH = [];
if ~isnan(x1d.b2OneXTransformFH)
    xT1FH = str2func(x1d.b2OneXTransformFH);
end
if ~isnan(x1d.b2TwoXTransformFH)
    xT2FH = str2func(x1d.b2TwoXTransformFH);
end

%% load data
fName = fullfile(rc.pathStr,sprintf(rc.dataPat, uo.Subject, uo.DateStr));
[bs ds] = getBehavDataForDay(...
    'Filename', fName, ...
    'DataIndex', x1d.DataBlock, ...
    'WhichTrials', whichTrs, ...
    'DoCorrectEarlies', x1d.DoCorrectEarlies, ...
    'B2OneXTransformFH', xT1FH, ...
    'B2TwoXTransformFH', xT2FH, ...
    'MergeBlock1And2', mergeB1And2);

if nargout > 6
    tMatFilename = fitX1d.MatFilename;
    
    % kludge to paper over bug
    [tPath tName tExt] = fileparts(tMatFilename);
    if strcmp(tExt, '.pdf')
        tMatFilename = fullfile(strrep(tPath, '/pdfFits', '/fitMatStats'), ...
            [tName, '.mat']);
    end

    load(tMatFilename, 'fitS', 'bootS');
end