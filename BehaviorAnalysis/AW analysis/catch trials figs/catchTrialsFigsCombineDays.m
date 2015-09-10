clear all;
rc.pathStr = 'Y:\home\andrew\Behavior\Data';
rc.dataPat = 'data-i%03d-%s.mat';
rootDir = 'C:\Users\ashley\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis\AW analysis\';
% rc.indexFilename = fullfile(rootDir, 'experimentIndexes\subj-days-lg.xls');
rc.indexcatchFilename = fullfile(rootDir, 'catch trials figs\catch-days.xlsx');
% rc.fitOutputFilename = fullfile(rootDir, 'experimentIndexes\subj-fits-lg.xls');
% rc.fitOutputSummary = fullfile(rootDir, 'output\summary');
% rc.fitOutputPdfDir = fullfile(rootDir, 'output\pdfFits');
% rc.fitOutputMatDir = fullfile(rootDir, 'output\fitMatStats');

%load xls file and data path
[Ldata, Ltext, Lraw] = xlsread(rc.indexcatchFilename);
cd(rc.pathStr);

%load selected trials from first data struct
ID = Ldata(1,1);
dayStr = num2str(Ldata(1,2));
daysAll = unique(Ldata(:,2));
dayStrRange = [num2str(daysAll(1,1)) '-' num2str(daysAll(end,1))];
time = num2str(Ldata(1,3));
if str2num(time) < 1000
    timeStr = num2str(time);
    timeStr = ['0', timeStr];
else
    timeStr = num2str(time);
end
trRange = Ltext(2,4); %if you do have a trial range
trStr = trRange{1};
dataPat = 'data-i%03d-%s-%s.mat';
fitdataName = sprintf(dataPat, ID, dayStr,timeStr);
load(fitdataName);



        
m = any(unique(cell2mat_padded(input.tCatchGratingDirectionDeg)));
    if m == 0
        %auditory 
        atrialOutcomeCell = eval(['input.trialOutcomeCell(', trStr,')']);
        atBlock2TrialNumber = eval(['input.tBlock2TrialNumber(', trStr,')']);
        atCatchGratingDirectionDeg = eval(['input.tCatchGratingDirectionDeg(', trStr,')']);
        atGratingDirectionDeg = eval(['input.tGratingDirectionDeg(', trStr,')']);
        atGratingContrast = eval(['input.tGratingContrast(', trStr,')']);
        acatchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')']);
        atCyclesOn = eval(['input.tCyclesOn(', trStr,')']);
        atSoundCatchAmplitude = eval(['input.tSoundCatchAmplitude(', trStr,')']);
        atSoundTargetAmplitude = eval(['input.tSoundTargetAmplitude(', trStr,')']);
        acatchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')']);
        trialOutcomeCell = [];
        tBlock2TrialNumber = [];
        tCatchGratingDirectionDeg = [];
        tGratingDirectionDeg = [];
        tGratingContrast = [];
        catchTrialOutcomeCell = [];
        tCyclesOn = [];
    else 
        %visual
        trialOutcomeCell = eval(['input.trialOutcomeCell(', trStr,')']);
        tBlock2TrialNumber = eval(['input.tBlock2TrialNumber(', trStr,')']);
        tCatchGratingDirectionDeg = eval(['input.tCatchGratingDirectionDeg(', trStr,')']);
        tGratingDirectionDeg = eval(['input.tGratingDirectionDeg(', trStr,')']);
        tGratingContrast = eval(['input.tGratingContrast(', trStr,')']);
        catchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')']);
        tCyclesOn = eval(['input.tCyclesOn(', trStr,')']);
        atrialOutcomeCell = [];
        atBlock2TrialNumber = [];
        atCatchGratingDirectionDeg = [];
        atGratingDirectionDeg = [];
        atGratingContrast = [];
        acatchTrialOutcomeCell = [];
        atCyclesOn = [];
        atSoundCatchAmplitude = [];
        atSoundTargetAmplitude = [];
        acatchTrialOutcomeCell = [];

    end
    
    if size(trialOutcomeCell,1) > size(trialOutcomeCell,2)
        atrialOutcomeCell = eval(['input.trialOutcomeCell(', trStr,')'])';
        atBlock2TrialNumber = eval(['input.tBlock2TrialNumber(', trStr,')'])';
        atCatchGratingDirectionDeg = eval(['input.tCatchGratingDirectionDeg(', trStr,')'])';
        atGratingDirectionDeg = eval(['input.tGratingDirectionDeg(', trStr,')'])';
        atGratingContrast = eval(['input.tGratingContrast(', trStr,')'])';
        acatchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')'])';
        atCyclesOn = eval(['input.tCyclesOn(', trStr,')'])';
        atSoundCatchAmplitude = eval(['input.tSoundCatchAmplitude(', trStr,')'])';
        atSoundTargetAmplitude = eval(['input.tSoundTargetAmplitude(', trStr,')'])';
        acatchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')'])';
        trialOutcomeCell = eval(['input.trialOutcomeCell(', trStr,')'])';
        tBlock2TrialNumber = eval(['input.tBlock2TrialNumber(', trStr,')'])';
        tCatchGratingDirectionDeg = eval(['input.tCatchGratingDirectionDeg(', trStr,')'])';
        tGratingDirectionDeg = eval(['input.tGratingDirectionDeg(', trStr,')'])';
        tGratingContrast = eval(['input.tGratingContrast(', trStr,')'])';
        catchTrialOutcomeCell = eval(['input.catchTrialOutcomeCell(', trStr,')'])';
        tCyclesOn = eval(['input.tCyclesOn(', trStr,')'])';
    end
        
    
    
    




%add trials from additional structs
daysDouble = size(Ldata(:,1));
nDays = daysDouble(1);
for i = 2:nDays
    ID = Ldata(i,1);
    dayStr = num2str(Ldata(i,2));
   
    time = num2str(Ldata(i,3));
    if str2num(time) < 1000
        timeStr = num2str(time);
        timeStr = ['0', timeStr];
    else
        timeStr = num2str(time);
    end
    trRange = Ltext(i+1,4);
    trStr = trRange{1};
    dataPat = 'data-i%03d-%s-%s.mat';
    fitdataName = sprintf(dataPat, ID, dayStr,timeStr);
    load(fitdataName);
    b = any(unique(cell2mat_padded(input.tCatchGratingDirectionDeg)));
    if b == 0
        try
        atrialOutcomeCell = [atrialOutcomeCell eval(['input.trialOutcomeCell(', trStr,')'])];
        atBlock2TrialNumber = [atBlock2TrialNumber eval(['input.tBlock2TrialNumber(', trStr,')'])];
        atCatchGratingDirectionDeg = [atCatchGratingDirectionDeg eval(['input.tCatchGratingDirectionDeg(', trStr,')'])];
        atGratingDirectionDeg = [atGratingDirectionDeg eval(['input.tGratingDirectionDeg(', trStr,')'])];
        atGratingContrast = [atGratingContrast eval(['input.tGratingContrast(', trStr,')'])];
        acatchTrialOutcomeCell = [acatchTrialOutcomeCell eval(['input.catchTrialOutcomeCell(', trStr,')'])];
        atCyclesOn = [atCyclesOn eval(['input.tCyclesOn(', trStr,')'])];
        atSoundCatchAmplitude = [atSoundCatchAmplitude eval(['input.tSoundCatchAmplitude(', trStr,')'])];
        atSoundTargetAmplitude = [atSoundTargetAmplitude eval(['input.tSoundTargetAmplitude(', trStr,')'])];
        catch
        atrialOutcomeCell = [atrialOutcomeCell eval(['input.trialOutcomeCell(', trStr,')'])'];
        atBlock2TrialNumber = [atBlock2TrialNumber eval(['input.tBlock2TrialNumber(', trStr,')'])'];
        atCatchGratingDirectionDeg = [atCatchGratingDirectionDeg eval(['input.tCatchGratingDirectionDeg(', trStr,')'])'];
        atGratingDirectionDeg = [atGratingDirectionDeg eval(['input.tGratingDirectionDeg(', trStr,')'])'];
        atGratingContrast = [atGratingContrast eval(['input.tGratingContrast(', trStr,')'])'];
        acatchTrialOutcomeCell = [acatchTrialOutcomeCell eval(['input.catchTrialOutcomeCell(', trStr,')'])'];
        atCyclesOn = [atCyclesOn eval(['input.tCyclesOn(', trStr,')'])'];
        atSoundCatchAmplitude = [atSoundCatchAmplitude eval(['input.tSoundCatchAmplitude(', trStr,')'])'];
        atSoundTargetAmplitude = [atSoundTargetAmplitude eval(['input.tSoundTargetAmplitude(', trStr,')'])'];
        end
        
    else
        try
        trialOutcomeCell = [trialOutcomeCell eval(['input.trialOutcomeCell(', trStr,')'])];
        tBlock2TrialNumber = [tBlock2TrialNumber eval(['input.tBlock2TrialNumber(', trStr,')'])];
        tCatchGratingDirectionDeg = [tCatchGratingDirectionDeg eval(['input.tCatchGratingDirectionDeg(', trStr,')'])];            
        tGratingDirectionDeg = [tGratingDirectionDeg eval(['input.tGratingDirectionDeg(', trStr,')'])];
        tGratingContrast = [tGratingContrast eval(['input.tGratingContrast(', trStr,')'])];
        catchTrialOutcomeCell = [catchTrialOutcomeCell eval(['input.catchTrialOutcomeCell(', trStr,')'])];
        tCyclesOn = [tCyclesOn eval(['input.tCyclesOn(', trStr,')'])];
        catch
        trialOutcomeCell = [trialOutcomeCell eval(['input.trialOutcomeCell(', trStr,')'])'];
        tBlock2TrialNumber = [tBlock2TrialNumber eval(['input.tBlock2TrialNumber(', trStr,')'])'];
        tCatchGratingDirectionDeg = [tCatchGratingDirectionDeg eval(['input.tCatchGratingDirectionDeg(', trStr,')'])'];            
        tGratingDirectionDeg = [tGratingDirectionDeg eval(['input.tGratingDirectionDeg(', trStr,')'])'];
        tGratingContrast = [tGratingContrast eval(['input.tGratingContrast(', trStr,')'])'];
        catchTrialOutcomeCell = [catchTrialOutcomeCell eval(['input.catchTrialOutcomeCell(', trStr,')'])'];
        tCyclesOn = [tCyclesOn eval(['input.tCyclesOn(', trStr,')'])'];            
        end
    end
end


%save location
fnout = ['Z:\Analysis\AW13\behavior\short catch experiments\summary figs'];
cd(fnout)

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);


%graph1 auditory catch
    Allcycles = cell2mat_padded(atCyclesOn);
    Cycles = unique(Allcycles);
    A = find(Allcycles == 1);
    atrialOutcomeCell(A)= [];
    atSoundCatchAmplitude(A) = [];
    atGratingDirectionDeg(A)=[];
    atSoundTargetAmplitude(A) = [];
    acatchTrialOutcomeCell(A) = [];
    atBlock2TrialNumber(A)= [];
    hitInd = strcmp(atrialOutcomeCell,'success');
    falsealarmInd = strcmp(atrialOutcomeCell,'failure');
    missInd = strcmp(atrialOutcomeCell,'ignore');
    block2 = cell2mat_padded(atBlock2TrialNumber);
    visInd = find(block2 == 1);
    audInd = find(block2 == 0);
    direction = celleqel2mat_padded(atSoundTargetAmplitude);
    dirs = unique(direction);
    doCatch = cell2mat_padded(atSoundCatchAmplitude);
    catchDeg = unique(doCatch);
    catchYesInd = strcmp(acatchTrialOutcomeCell,'FA');
    catchFalseAlarmInd = strcmp(acatchTrialOutcomeCell,'failure');
    catchTrial = doCatch > 0;
    
    for i = 1:length(dirs)-1
        trials = find(direction == dirs(i+1));
        visHitRate(i) = sum(hitInd(trials))/(sum(hitInd(trials))+sum(missInd(trials)));
        trialN(i) = sum(hitInd(trials))+sum(missInd(trials));
    end

    for i = 1:length(catchDeg)-1
        trials = find(doCatch == catchDeg(i+1));
        catchHitRate(i) = sum(catchYesInd(trials))/(sum(catchTrial(trials)));
        catchTrialN(i) = sum(catchTrial(trials));
    end
    

    %total number of trials
    nVisTrials = sum(hitInd(block2 ==0))+sum(missInd)
    nCatchTrials = sum(catchTrial)



    for i = 1:length(dirs)-1
        for iboot = 1:1000
            trials = find(direction == dirs(i+1));
            trialIndex = randi([1 length(trials)],length(trials),1);
            bootstrapIndex = trials(trialIndex);
            visHitRateBoots(i,iboot) = sum(hitInd(bootstrapIndex))/(sum(hitInd(bootstrapIndex))+sum(missInd(bootstrapIndex)));
        end  
    end
    visHitRateCI = prctile(visHitRateBoots,[2.5 97.5],2);

    for i = 1:length(catchDeg)-1
        for iboot = 1:1000
            trials = find(doCatch == catchDeg(i+1));
            trialIndex = randi([1 length(trials)],length(trials),1);
            bootstrapIndex = trials(trialIndex);
        catchHitRateBoots(i,iboot) = sum(catchYesInd(bootstrapIndex))/(sum(catchTrial(bootstrapIndex)));
        end
    end
    catchHitRateCI = prctile(catchHitRateBoots,[2.5 97.5],2);
    
    xLimL = [min(dirs(~dirs==0))*0.5, max(dirs)*1.5];
    if length(xLimL)==1, 
        xLimL = get(gca, 'XLim'); 
    end
    xL1 = [floor(log10(xLimL(1))) ceil(log10(xLimL(2)))];
    xTickL = 10.^(xL1(1):1:xL1(2));
    xTickL = xTickL(xTickL>=xLimL(1) & xTickL<=xLimL(2));
    if length(xTickL) == 0, 
        % no log10 ticks in range, create two and set them
        xTickL = [xLimL(1) xLimL(2)]; %[xTL(1)/10 xTL(1) xTL(1)]
        xTickL = chop(xTickL,2);
    end
        
        xTickL = sort([xTickL dirs]);
        xTLabelL = cellstr(num2str(xTickL(:)));
    if xLimL(1) > -2
        xLimL(1) = -2;
    end
    
    figure;
    errorbar(dirs(2:end),visHitRate,visHitRate-visHitRateCI(:,1)',visHitRateCI(:,2)'-visHitRate,'go','LineStyle', 'none');
    trialNLabel = strcat(repmat({'\leftarrow'}, [length(trialN) 1]),cellstr(num2str(trialN(:))), repmat({' trials'}, [length(trialN) 1]))';
    text(dirs(2:end),visHitRate,trialNLabel)
    % errorbar(dirs(2:end),visHitRate,visHitRate-visHitRateCI(:,1)',visHitRateCI(:,2)'-visHitRate,'g','Linewidth',3);
    hold on
    % plot(0,audHitRate,'.-','Linewidth',5);
    % hold on
    errorbar(catchDeg(2:end),catchHitRate,catchHitRate-catchHitRateCI(:,1)',catchHitRateCI(:,2)'-catchHitRate,'ko','LineStyle','none');
    catchTrialNLabel = strcat(repmat({'\leftarrow'}, [length(catchTrialN) 1]),cellstr(num2str(catchTrialN(:))), repmat({' trials'}, [length(catchTrialN) 1]))';
    text(catchDeg(2:end),catchHitRate,catchTrialNLabel)
    % errorbar(catchDeg(2:end),catchHitRate,catchHitRate-catchHitRateCI(:,1)',catchHitRateCI(:,2)'-catchHitRate,'k','Linewidth',3);
    hold on
%     xlim([0 0.25]);
%     ylim([1 100]);
    set(gca,'xscale','log');
    ylabel('Hits/Hits+Miss')
    xlabel('Sound Amplitude')
    title([num2str(ID) ' - Aud Catch Trials ' dayStrRange]);
%     legendInfo = {num2str(nVisTrials), num2str(nCatchTrials)};
%     legend(legendInfo)
    hold on
    % set(gca, 'XGrid', 'on');
    % set(gca, 'XTick', xTickL); 
    set(gca, 'XTick', [1 10 100]); 
    % set(gca, 'XTickLabel', xTLabelL);
    % for i = catchDeg
    %     vline(i,'k');
    %     hold on
    % end

print([fnout ['\' [num2str(ID) ' - Aud Catch Trials ' dayStrRange] '.pdf']], '-dpdf')    
    
    %graph2
    Allcycles = cell2mat_padded(tCyclesOn);
    Cycles = unique(Allcycles);
    A = find(Allcycles == 1);
    trialOutcomeCell(A)= [];
    tCatchGratingDirectionDeg(A) = [];
    tGratingDirectionDeg(A)=[];
    catchTrialOutcomeCell(A) = [];
    tBlock2TrialNumber(A)= [];
    hitInd = strcmp(trialOutcomeCell,'success');
    falsealarmInd = strcmp(trialOutcomeCell,'failure');
    missInd = strcmp(trialOutcomeCell,'ignore');
    block2 = cell2mat_padded(tBlock2TrialNumber);
    visInd = find(block2 == 0);
    audInd = find(block2 == 1);
    direction = cell2mat_padded(tGratingDirectionDeg);
    dirs = unique(direction);
    doCatch = cell2mat_padded(tCatchGratingDirectionDeg);
    catchDeg = unique(doCatch);
    catchYesInd = strcmp(catchTrialOutcomeCell,'FA');
    catchFalseAlarmInd = strcmp(catchTrialOutcomeCell,'failure');
    catchTrial = doCatch > 0;


    for i = 1:length(dirs)-1
        trials = find(direction == dirs(i+1));
        visHitRate1(i) = sum(hitInd(trials))/(sum(hitInd(trials))+sum(missInd(trials)));
        trialN1(i) = sum(hitInd(trials))+sum(missInd(trials));
    end

    for i = 1:length(catchDeg)-1
        trials = find(doCatch == catchDeg(i+1));
        catchHitRate1(i) = sum(catchYesInd(trials))/(sum(catchTrial(trials)));
        catchTrialN1(i) = sum(catchTrial(trials));
    end

    %total number of trials
    nVisTrials = sum(hitInd(block2 ==0))+sum(missInd)
    nCatchTrials = sum(catchTrial)


visHitRateBoots = zeros(length(dirs)-1,1000);
    for i = 1:length(dirs)-1
        for iboot = 1:1000
            trials = find(direction == dirs(i+1));
            trialIndex = randi([1 length(trials)],length(trials),1);
            bootstrapIndex = trials(trialIndex);
            visHitRateBoots(i,iboot) = sum(hitInd(bootstrapIndex))/(sum(hitInd(bootstrapIndex))+sum(missInd(bootstrapIndex)));
        end  
    end
    visHitRateCI1 = prctile(visHitRateBoots,[2.5 97.5],2);
    
    
    catchHitRateBoots = zeros(length(catchDeg)-1,1000);
    for i = 1:length(catchDeg)-1
        for iboot = 1:1000
            trials = find(doCatch == catchDeg(i+1));
            trialIndex = randi([1 length(trials)],length(trials),1);
            bootstrapIndex = trials(trialIndex);
        catchHitRateBoots(i,iboot) = sum(catchYesInd(bootstrapIndex))/(sum(catchTrial(bootstrapIndex)));
        end
    end
    catchHitRateCI1 = prctile(catchHitRateBoots,[2.5 97.5],2);
    
    
xLimL = [min(dirs(~dirs==0))*0.5, max(dirs)*1.5];
if length(xLimL)==1, 
    xLimL = get(gca, 'XLim'); 
end
xL1 = [floor(log10(xLimL(1))) ceil(log10(xLimL(2)))];
xTickL = 10.^(xL1(1):1:xL1(2));
xTickL = xTickL(xTickL>=xLimL(1) & xTickL<=xLimL(2));
if length(xTickL) == 0, 
    % no log10 ticks in range, create two and set them
    xTickL = [xLimL(1) xLimL(2)]; %[xTL(1)/10 xTL(1) xTL(1)]
    xTickL = chop(xTickL,2);
end
    xTickL = sort([xTickL dirs']);
    xTLabelL = cellstr(num2str(xTickL(:)));
if xLimL(1) > -2
    xLimL(1) = -2;
end
    
figure;
errorbar(dirs(2:end),visHitRate1,visHitRate1-visHitRateCI1(:,1)',visHitRateCI1(:,2)'-visHitRate1,'go','LineStyle', 'none');
trialNLabel = strcat(repmat({'\leftarrow'}, [length(trialN1) 1]),cellstr(num2str(trialN1(:))), repmat({' trials'}, [length(trialN1) 1]))';
text(dirs(2:end),visHitRate1,trialNLabel)
% errorbar(dirs(2:end),visHitRate,visHitRate-visHitRateCI(:,1)',visHitRateCI(:,2)'-visHitRate,'g','Linewidth',3);
hold on
% plot(0,audHitRate,'.-','Linewidth',5);
% hold on
errorbar(catchDeg(2:end),catchHitRate1,catchHitRate1-catchHitRateCI1(:,1)',catchHitRateCI1(:,2)'-catchHitRate1,'ko','LineStyle','none');
catchTrialNLabel = strcat(repmat({'\leftarrow'}, [length(catchTrialN1) 1]),cellstr(num2str(catchTrialN1(:))), repmat({' trials'}, [length(catchTrialN1) 1]))';
text(catchDeg(2:end),catchHitRate1,catchTrialNLabel)
% errorbar(catchDeg(2:end),catchHitRate,catchHitRate-catchHitRateCI(:,1)',catchHitRateCI(:,2)'-catchHitRate,'k','Linewidth',3);
hold on
ylim([0 1]);
xlim([1 100]);
set(gca,'xscale','log');
ylabel('Hits/Hits+Miss')
xlabel('Deg Ori Change')
title([num2str(ID) ' - Vis Catch Trials ' dayStrRange]);
%     legendInfo = {num2str(nVisTrials), num2str(nCatchTrials)};
%     legend(legendInfo)
hold on
% set(gca, 'XGrid', 'on');
% set(gca, 'XTick', xTickL); 
set(gca, 'XTick', [1 10 100]); 
% set(gca, 'XTickLabel', xTLabelL);
% for i = catchDeg
%     vline(i,'k');
%     hold on
% end

print([fnout ['\' [num2str(ID) ' - Vis Catch Trials ' dayStrRange] '.pdf']], '-dpdf')    




