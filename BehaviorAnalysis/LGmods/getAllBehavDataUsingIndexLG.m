function [bs ds x1d fitX1d xd fitXd fitS bootS] = getAllBehavDataUsingIndexLG(varargin)
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
    'DataBlock', 'max', ...
    'LoadDataDebug', false, ...
    'FitXd', []};
    
uo = stropt2struct(stropt_defaults(userDefs, varargin));

%% get xls data and find params
rc = behavConstsHADC8;
xd = frm_xls2frm(rc.indexFilename);

% do explicit casting for cases where xls may have different types
if all(isnumeric(xd.DateStr))
    xd.DateStr = cellstr(int2str(xd.DateStr(:)))';
end

if isempty(uo.FitXd)
    fitXd = frm_xls2frm(rc.fitOutputFilename);
    if all(isnumeric(fitXd.DateStr))
        fitXd.DateStr = cellstr(int2str(fitXd.DateStr(:)))';
    end
    if isnumeric(fitXd.DateTimeStarted)
        assert(all(isnan(fitXd.DateTimeStarted)), ...
            'if numeric must be all empty - timestamp is str');
        fitXd.DateTimeStarted = cell([1 fitXd.nRows]);
    end
else
    fitXd = uo.FitXd;
end

%% find index row
xdN = frm_findrownum(xd, ...
    { 'Subject', uo.Subject, ...
    'DateStr', uo.DateStr, ...
    'DataBlock', uo.DataBlock, ...
    });
x1d = frm_extractrow(xd, xdN);
if isnan(x1d.DataBlock), x1d.DataBlock = 'max'; end

% fit row
fitRowN = frm_findrownum(fitXd, ...
    { 'Subject', uo.Subject, ...
    'DateStr', uo.DateStr, ...
    'DataBlock', uo.DataBlock }, ...
    'IgnoreMissing', true);
if isnan(fitRowN) 
    fitX1d = [];
else
    fitX1d = frm_extractrow(fitXd, fitRowN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse whichtrials
whichTrs = x1d.TrialRangeToUse;
if isnan(whichTrs)
    whichTrs = [];
end
if length(whichTrs)== 3
    if whichTrs == 'NaN';
        whichTrs = [];
    end
end
if ~isempty(whichTrs)
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
if isempty(x1d.DoBlock1Only)
    doB1only = false;
else
    doB1only = x1d.DoBlock1Only;
end
if isempty(x1d.SplitBlock1)
    splitB1 = false;
else
    splitB1 = x1d.SplitBlock1;
end
if isempty(x1d.MergeAllMatFiles)
    mergeMats = false;
else
    mergeMats = x1d.MergeAllMatFiles;
end
%% load data
fName = fullfile(rc.pathStr,sprintf(rc.dataPat, uo.Subject, uo.DateStr));
[bs ds] = getBehavDataForDayLG(...
    'Filename', fName, ...
    'DataIndex', uo.DataBlock, ...
    'WhichTrials', whichTrs, ...
    'DoCorrectEarlies', x1d.DoCorrectEarlies, ...
    'MergeBlock1And2', mergeB1And2, ...
    'DoBlock1Only', doB1only, ...
    'SplitBlock1', splitB1, ...
    'MergeMats', mergeMats);

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