%% vis catch trials

subnum = 'i613';
dateofexp = '151015';
timeofexp = '1114';

load(fullfile('Y:\home\andrew\Behavior\Data',['data-' subnum '-' dateofexp '-' timeofexp]))

% get catch trial hit rates
catchGrating = celleqel2mat_padded(input.tCatchGratingDirectionDeg);
catchOutcome = input.catchTrialOutcomeCell;
cDirs = unique(catchGrating(~isnan(catchGrating)));
cind = [];
cN = zeros(1,length(cDirs));
faN = zeros(1,length(cDirs));
crN = zeros(1,length(cDirs));
for idir = 1:length(cDirs)
    cind{idir} = find(catchGrating == cDirs(idir));
    cN(idir) = length(find(catchGrating == cDirs(idir)));
    faN(idir) = length(intersect(find(strcmp(catchOutcome,'FA')),(cind{idir})));
    criN(idir) = length(intersect(find(strcmp(catchOutcome,'CR')),(cind{idir})));
end

cHR_vis = faN./cN;

%plot on same axis as valid trial data
grating = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(grating(~isnan(grating)));

minI = dirs(2);
maxI = dirs(end);

xgrid = logspace(log10(minI*0.5),log10(maxI*1.5),100);
xLim = [min(xgrid) max(xgrid)] .* 1 .* [0.75 1.25];
xLim = 10.^ceil(log10(xLim) - [1 0]);

exampleCatchFigs = figure;
subplot(2,1,1)
plot(cDirs,cHR_vis,'ko-')

set(gca, 'XLim', xLim);
set(gca, 'XScale', 'log');
set(gca, 'YLim', [0 1]);

xTick = get(gca, 'XTick');
xTL = eval(['{' sprintf('''%g'',', chop(xTick,2)) '}']);
set(gca, 'XTickLabel', xTL);
xlabel('orientation');
ylabel('fraction correct');
title([num2str(subnum) ' ' num2str(dateofexp)])


%% aud catch trials

subnum = 'i613';
dateofexp = '150922';
timeofexp = '1002';

load(fullfile('Y:\home\andrew\Behavior\Data',['data-' subnum '-' dateofexp '-' timeofexp]))

% get catch trial hit rates
catchAmplitude = celleqel2mat_padded(input.tSoundCatchAmplitude);
catchOutcome = input.catchTrialOutcomeCell;
cAmps = unique(catchAmplitude(~isnan(catchAmplitude)));
cind = [];
cN = zeros(1,length(cAmps));
faN = zeros(1,length(cAmps));
crN = zeros(1,length(cAmps));
for iamp = 1:length(cAmps)
    cind{iamp} = find(catchAmplitude == cAmps(iamp));
    cN(iamp) = length(find(catchAmplitude == cAmps(iamp)));
    faN(iamp) = length(intersect(find(strcmp(catchOutcome,'FA')),(cind{iamp})));
    criN(iamp) = length(intersect(find(strcmp(catchOutcome,'CR')),(cind{iamp})));
end

cHR_aud = faN./cN;

%plot on same axis as valid trial data
amplitude = celleqel2mat_padded(input.tSoundTargetAmplitude);
amps = unique(amplitude(~isnan(amplitude)));

minI = amps(2);
maxI = amps(end);

xgrid = logspace(log10(minI*0.5),log10(maxI*1.5),100);
xLim = [min(xgrid) max(xgrid)] .* 1 .* [0.75 1.25];
xLim = 10.^ceil(log10(xLim) - [1 0]);

figure(exampleCatchFigs);
subplot(2,1,2)
plot(cAmps,cHR_aud,'ko-')

set(gca, 'XLim', xLim);
set(gca, 'XScale', 'log');
set(gca, 'YLim', [0 1]);

xTick = get(gca, 'XTick');
xTL = eval(['{' sprintf('''%g'',', chop(xTick,2)) '}']);
set(gca, 'XTickLabel', xTL);
xlabel('volume');
ylabel('fraction correct');
title([num2str(subnum) ' ' num2str(dateofexp)])

suptitle('catch trials')
%% save pdf
outName = fullfile('Z:\Analysis\Behavior\output\pdfFits',[subnum '-' dateofexp '-invalTrials'])
exportfig_print(gcf, outName, ...
        'FileFormat', 'pdf', ...
        'Size', [6 4.8+4]);




