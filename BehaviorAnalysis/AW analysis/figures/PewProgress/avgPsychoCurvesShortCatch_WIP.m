%% data 1
subNum = '613';
date = '150403';
time = '1046';

CD = 'Y:\home\andrew\Behavior\Data';
cd(CD);

mworks = ['data-i' subNum '-' date '-' time '.mat'];
load(mworks)

hitInd = strcmp(input.trialOutcomeCell,'success');
falsealarmInd = strcmp(input.trialOutcomeCell,'failure');
missInd = strcmp(input.trialOutcomeCell,'ignore');
block2 = cell2mat_padded(input.tBlock2TrialNumber);
visInd = find(block2 == 0);
audInd = find(block2 == 1);
direction = cell2mat_padded(input.tGratingDirectionDeg);
dirs = unique(direction);
doCatch = cell2mat_padded(input.tCatchGratingDirectionDeg);
catchDeg = unique(doCatch);
catchYesInd = strcmp(input.catchTrialOutcomeCell,'FA');
catchFalseAlarmInd = strcmp(input.catchTrialOutcomeCell,'failure');
catchTrial = doCatch > 0;

%auditory trials


clear input
%% data 2
date = '150318';
time = '0926';

CD = 'Y:\home\andrew\Behavior\Data';
cd(CD);

mworks = ['data-i' subNum '-' date '-' time '.mat'];
load(mworks)

hitInd = cat(2,hitInd,strcmp(input.trialOutcomeCell,'success'));
falsealarmInd = cat(2,falsealarmInd,strcmp(input.trialOutcomeCell,'failure'));
missInd = cat(2,missInd,strcmp(input.trialOutcomeCell,'ignore'));
block2 = cat(1,block2,cell2mat_padded(input.tBlock2TrialNumber));
% visInd = find(block2 == 0);
% audInd = find(block2 == 1);
direction = cat(1,direction,cell2mat_padded(input.tGratingDirectionDeg));
% dirs = unique(direction);
doCatch = cat(1,doCatch,cell2mat_padded(input.tCatchGratingDirectionDeg));
% catchDeg = unique(doCatch);
catchYesInd = cat(2,catchYesInd,strcmp(input.catchTrialOutcomeCell,'FA'));
% catchFalseAlarmInd = strcmp(input.catchTrialOutcomeCell,'failure');
catchTrial = cat(1,catchTrial,(cell2mat_padded(input.tCatchGratingDirectionDeg)) > 0);

%% data 3 - cropped
date = '150326';
time = '1106';

CD = 'Y:\home\andrew\Behavior\Data';
cd(CD);

mworks = ['data-i' subNum '-' date '-' time '.mat'];
load(mworks)

hitIndcrop = strcmp(input.trialOutcomeCell,'success');
hitIndcrop = hitIndcrop(101:500);
falsealarmIndcrop = strcmp(input.trialOutcomeCell,'failure');
falsealarmIndcrop = falsealarmIndcrop(101:500);
missIndcrop = strcmp(input.trialOutcomeCell,'ignore');
missIndcrop = missIndcrop(101:500);
block2crop = cell2mat_padded(input.tBlock2TrialNumber);
block2crop = block2crop(101:500);
directioncrop = cell2mat_padded(input.tGratingDirectionDeg);
directioncrop = directioncrop(101:500);
doCatchcrop = cell2mat_padded(input.tCatchGratingDirectionDeg);
doCatchcrop = doCatchcrop(101:500);
catchYesIndcrop = strcmp(input.catchTrialOutcomeCell,'FA');
catchYesIndcrop = catchYesIndcrop(101:500);
catchFalseAlarmIndcrop = strcmp(input.catchTrialOutcomeCell,'failure');
catchFalseAlarmIndcrop = catchFalseAlarmIndcrop(101:500);
catchTrialcrop = cell2mat_padded(input.tCatchGratingDirectionDeg) > 0;
catchTrialcrop = catchTrialcrop(101:500);

hitInd = cat(2,hitInd,hitIndcrop);
falsealarmInd = cat(2,falsealarmInd,falsealarmIndcrop);
missInd = cat(2,missInd,missIndcrop);
block2 = cat(1,block2,block2crop);
visInd = find(block2 == 0);
audInd = find(block2 == 1);
direction = cat(1,direction,directioncrop);
dirs = unique(direction);
doCatch = cat(1,doCatch,doCatchcrop);
catchDeg = unique(doCatch);
catchYesInd = cat(2,catchYesInd,catchYesIndcrop);
catchTrial = cat(1,catchTrial,catchTrialcrop);

%% calc
for i = 1:length(dirs)-1
    trials = find(direction == dirs(i+1));
    visHitRate(i) = sum(hitInd(trials))/(sum(hitInd(trials))+sum(missInd(trials)));
    [binoVisHitRate(i) binoVis(i,:)] = binofit(sum(hitInd(trials)),(sum(hitInd(trials))+sum(missInd(trials))));
end

for i = 1:length(catchDeg)-1
    trials = find(doCatch == catchDeg(i+1));
    catchHitRate(i) = sum(catchYesInd(trials))/(sum(catchTrial(trials)));
    [binoCatchHitRate(i) binoCatch(i,:)] = binofit(sum(catchYesInd(trials)),(sum(catchTrial(trials))));
end



%auditory trials
trials = find(block2 == 1);
audHitRate = sum(hitInd(trials))/(sum(hitInd(trials))+sum(missInd(trials)));
[binoAudHitRate binoAud] = binofit(sum(hitInd(trials)),(sum(hitInd(trials))+sum(missInd(trials))));


%total number of trials
nVisTrials = sum(hitInd(block2 ==0))+sum(missInd)
nCatchTrials = sum(catchTrial)
nAudTrials = sum(hitInd(block2 == 1))
%%
binoVis = binofit(visHitRate,

%% calc bootstraps
for i = 1:length(dirs)-1
    for iboot = 1:1000
        trials = find(direction == dirs(i+1));
        trialIndex = randi([1 length(trials)],length(trials),1);
        bootstrapIndex = trials(trialIndex);
        visHitRateBoots(i,iboot) = sum(hitInd(bootstrapIndex))/(sum(hitInd(bootstrapIndex))+sum(missInd(bootstrapIndex)));
    end  
end
visHitRateCI = prctile(visHitRateBoots,[5 95],2);

for i = 1:length(catchDeg)-1
    for iboot = 1:1000
         trials = find(doCatch == catchDeg(i+1));
         trialIndex = randi([1 length(trials)],length(trials),1);
         bootstrapIndex = trials(trialIndex);
        catchHitRateBoots(i,iboot) = sum(catchYesInd(bootstrapIndex))/(sum(catchTrial(bootstrapIndex)));
    end
end
catchHitRateCI = prctile(catchHitRateBoots,[5 95],2);

for iboot = 1:1000
    trials = find(block2 == 1);
    trialIndex = randi([1 length(trials)],length(trials),1);
    bootstrapIndex = trials(trialIndex);
    audHitRateBoots(1,iboot) = sum(hitInd(bootstrapIndex))/(sum(hitInd(bootstrapIndex))+sum(missInd(bootstrapIndex)));
end
audHitRateCI = prctile(audHitRateBoots,[5 95],2);


%% PLOT
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
errorbar(dirs(2:end),visHitRate,visHitRate-visHitRateCI(:,1)',visHitRateCI(:,2)'-visHitRate,'go','LineStyle', 'none');
% errorbar(dirs(2:end),visHitRate,visHitRate-visHitRateCI(:,1)',visHitRateCI(:,2)'-visHitRate,'g','Linewidth',3);
hold on
% plot(0,audHitRate,'.-','Linewidth',5);
% hold on
errorbar(catchDeg(2:end),catchHitRate,catchHitRate-catchHitRateCI(:,1)',catchHitRateCI(:,2)'-catchHitRate,'ko','LineStyle','none');
% hold on
errorbar(100,audHitRate,audHitRate-audHitRateCI(:,1)',audHitRateCI(:,2)'-audHitRate,'ko','LineStyle','none');
% errorbar(catchDeg(2:end),catchHitRate,catchHitRate-catchHitRateCI(:,1)',catchHitRateCI(:,2)'-catchHitRate,'k','Linewidth',3);
hold on
ylim([0 1]);
xlim([1 100]);
set(gca,'xscale','log');
ylabel('Hits/Hits+Miss')
xlabel('Deg Ori Change')
title(['AW13 - Catch Trials' num2str(date)]);
hold on
% set(gca, 'XGrid', 'on');
% set(gca, 'XTick', xTickL); 
set(gca, 'XTick', [1 10 100]); 
% set(gca, 'XTickLabel', xTLabelL);
% for i = catchDeg
%     vline(i,'k');
%     hold on
% end

%% PLOT with binomial confidence intervals
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
errorbar(dirs(2:end),visHitRate,visHitRate-binoVis(:,1)',binoVis(:,2)'-visHitRate,'go','LineStyle', 'none');

hold on
errorbar(catchDeg(2:end),catchHitRate,catchHitRate-binoCatch(:,1)',binoCatch(:,2)'-catchHitRate,'ko','LineStyle','none');
hold on
errorbar(100,audHitRate,audHitRate-binoAud(:,1)',binoAud(:,2)'-audHitRate,'ko','LineStyle','none');
hold on
ylim([0 1]);
xlim([1 100]);
set(gca,'xscale','log');
ylabel('Hits/Hits+Miss')
xlabel('Deg Ori Change')
title(['AW13 - Catch Trials' num2str(date)]);
hold on
% set(gca, 'XGrid', 'on');
% set(gca, 'XTick', xTickL); 
set(gca, 'XTick', [1 10 100]); 
% set(gca, 'XTickLabel', xTLabelL);
% for i = catchDeg
%     vline(i,'k');
%     hold on
% end