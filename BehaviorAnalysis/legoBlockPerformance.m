function [ output_args ] = legoBlockPerformance(subjNum, date)
%LEGOBLOCKPERFORMANCE:: function to analyze average performance for the
%nth trial in a block of left or right trials for the Lego task. Must enter
%in subjNum and date in strings, with date being 'YYMMDD' format, with
%optional 'HHMM' suffix.

pathName = '~/Documents/MWorks/Data';
dataName = strcat('data-', subjNum, '-', date, '.mat');
fName = fullfile(pathName, dataName);

if ~exist(fName)
    n  = dir([fName(1:40) '*']);
    if size(n,1) == 1
        fName = fullfile(pathName, n.name);
        ds =  load(fName);
    elseif size(n,1) > 1
        error('Too many MAT files!!!');
    end
end
if isfield(ds, 'input')
    ds = ds.input;
end

%% Indexing and pre-analysis
tLeftTrial = cell2mat_padded(ds.tLeftTrial);
leftBlockN = cell2mat_padded(ds.tNBlockLeftTrsCompleted);
rightBlockN = cell2mat_padded(ds.tNBlockRightTrsCompleted);

successIx = strcmp(ds.trialOutcomeCell, 'success');
incorrectIx = strcmp(ds.trialOutcomeCell, 'incorrect');
ignoreIx = strcmp(ds.trialOutcomeCell, 'ignore');
as = struct;
subplotSize = {6,1};

%% Left Analysis
LBlockValues = unique(leftBlockN);
for i = 1:length(LBlockValues),
    trIx = (leftBlockN==LBlockValues(i)) & tLeftTrial;
    nTrs = sum(trIx);
    as.LnTrs(i) = nTrs;
    as.leftSuccess(i) = sum(trIx & successIx)/nTrs; 
    as.leftIncorrs(i) = sum(trIx & incorrectIx)/nTrs; 
    as.leftIgnore(i) = sum(trIx & ignoreIx)/nTrs; 
end
hold on
ax = subplot(subplotSize{:},1:2);
plot(LBlockValues, as.leftSuccess, 'k', LBlockValues, as.leftIncorrs, 'm', LBlockValues, as.leftIgnore, 'g')
xlabel('Block Position')
ylabel('Fraction of Outcomes')
title('Left Block Analysis')
legend('Success','Incorrect', 'Ignore', 'Location', 'Best')
xlim([1 length(LBlockValues)])
ax = subplot(subplotSize{:},3);
bar(LBlockValues, as.LnTrs)
ylabel('Trials per Block Position')
xlim([1 length(LBlockValues)])
%% Right Analysis
RBlockValues = unique(rightBlockN);
for i = 1:length(RBlockValues),
    trIx = (rightBlockN==RBlockValues(i)) & ~tLeftTrial;
    nTrs = sum(trIx);
    as.RnTrs(i) = nTrs;
    as.rightSuccess(i) = sum(trIx & successIx)/nTrs; 
    as.rightIncorrs(i) = sum(trIx & incorrectIx)/nTrs; 
    as.rightIgnore(i) = sum(trIx & ignoreIx)/nTrs; 
end
hold on
ax = subplot(subplotSize{:},4:5);
plot(RBlockValues, as.rightSuccess, 'k', RBlockValues, as.rightIncorrs, 'm', RBlockValues, as.rightIgnore, 'g')
xlabel('Block Position')
ylabel('Fraction of Outcomes')
title('Right Block Analysis')
legend('Success','Incorrect', 'Ignore', 'Location', 'Best')
xlim([1 length(RBlockValues)])
ax = subplot(subplotSize{:},6);
bar(RBlockValues, as.RnTrs)
ylabel('Trials per Block Position')
xlim([1 length(RBlockValues)])


%% Exporting Figure
figName = strcat('blockAnalysis-', subjNum, '-', date);
outputName = fullfile('~/Documents/MWorks/BehavOutputPdfs', figName);
exportfig_print(gcf, outputName , ...
             'FileFormat', 'pdf', ...
             'Size', [8 12], ...
             'PrintUI', false);
end