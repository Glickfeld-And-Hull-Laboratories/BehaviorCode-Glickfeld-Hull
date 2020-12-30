
SIx = s_cropped_include.SIx;
for i = 2:length(s_cropped_include); 
    success = s_cropped_include(i).SIx; 
    SIx = horzcat(SIx, success);
end 

FIx = s_cropped_include.FIx;
for i = 2:length(s_cropped_include); 
    failure = s_cropped_include(i).FIx; 
    FIx = horzcat(FIx, failure); 
end 

tLeftTrial = s_cropped_include.tLeftTrial;
for i = 2:length(s_cropped_include);
    TLT = s_cropped_include(i).tLeftTrial;
    tLeftTrial = horzcat(tLeftTrial, TLT); 
end

tGratingContrast = s_cropped_include.tGratingContrast;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).tGratingContrast;
    tGratingContrast = horzcat(tGratingContrast, TGC); 
end 

dGratingContrast = s_cropped_include.dGratingContrast;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).dGratingContrast;
    dGratingContrast = horzcat(dGratingContrast, TGC); 
end 

tGratingDiameter = s_cropped_include.tGratingDiameter;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).tGratingDiameter;
    tGratingDiameter = horzcat(tGratingDiameter, TGC); 
end 

dGratingDiameter = s_cropped_include.dGratingDiameter;
for i = 2:length(s_cropped_include); 
    DGC = s_cropped_include(i).dGratingDiameter;
    dGratingDiameter = horzcat(dGratingDiameter, DGC); 
end 

tRightResponse = s_cropped_include.tRR;
for i = 2:length(s_cropped_include); 
    tR_resp = s_cropped_include(i).tRR;
    tRightResponse = horzcat(tRightResponse, tR_resp); 
end 

tIgnore = s_cropped_include.tIgnore;
for i = 2:length(s_cropped_include); 
    Ignore = s_cropped_include(i).tIgnore;
    tIgnore = horzcat(tIgnore, Ignore); 
end 

tLeftResponse = zeros(length(tRightResponse),1) ;
tLeftResponse(intersect((find(tRightResponse == 0)),find(~tIgnore == 0))) = 1;

tBlock2 = s_cropped_include.tBlock2; 
for i = 2:length(s_cropped_include); 
    TBT = s_cropped_include(i).tBlock2; 
    tBlock2 = horzcat(tBlock2, TBT); 
end 

tDecisionTime = s_cropped_include.tDecisionTimeMs;
for i = 2:length(s_cropped_include); 
    TDT = s_cropped_include(i).tDecisionTimeMs;
    tDecisionTime = horzcat(tDecisionTime, TDT); 
end 

nTrials = length(tLeftTrial); 
tIgnore = zeros(1,nTrials); 
tIgnore(intersect(find(SIx == 0 ),find(FIx == 0))) = 1;
igno = find(tIgnore == 0) ;
% Rohan added 2 lines below
tGratingDiameterDeg=tGratingDiameter;
dGratingDiameterDeg=dGratingDiameter;
tGratingDiameter(find(tLeftTrial)) = tGratingDiameter(find(tLeftTrial)).*-1;

nTrials = length(tLeftTrial);
tLeftChoice = zeros(1,nTrials);
tLeftChoice(intersect(find(SIx),find(tLeftTrial))) = 1;
tLeftChoice(intersect(find(FIx),find(~tLeftTrial))) = 1;
%tLeftChoice now has successful left trials and unsuccessful not left
%trials, so shows all left responses both correct and incorrect 

SIx_sign = SIx;
FIx_sign = FIx;
SIx_sign(find(tLeftChoice)) = SIx(find(tLeftChoice)).*-1;
%this changes successful left trials to -1. Does not affect unsuccessful
%left trials because they are already zero. successful right trials = 1.

FIx_sign(find(tLeftChoice)) = FIx(find(tLeftChoice)).*-1;
%This changes unsuccessful left trials to -1. Does not affect successful
%left trials because they are alreay zero. unsuccessful right trials = 1.

tLeftChoice_prev = [NaN tLeftChoice(1:nTrials-1)];
SIx_prev = [NaN SIx_sign(1:nTrials-1)];
FIx_prev = [NaN FIx_sign(1:nTrials-1)];
% Assigns 1 as NAN (doesnt exist) and then assigns "previous trial values"
% for the 303 trials.



SIx = SIx(igno)';
FIx = FIx(igno)';
tGratingDiameterDeg = tGratingDiameterDeg(igno)'; 
dGratingDiameterDeg = dGratingDiameterDeg(igno)';
tBlock2 = tBlock2(igno)';
tIgnore = tIgnore(igno)';
tLeftTrial = tLeftTrial(igno)';
tRightResponse = tRightResponse(igno)';
tGratingDiameter = tGratingDiameter(igno)'; 
dGratingDiameter = dGratingDiameter(igno)';

tLeftTrial(find(tLeftTrial)) = -1;
tLeftTrial(find(~tLeftTrial)) = 1;

GSAA = abs(tGratingDiameterDeg)./dGratingDiameterDeg;
uGSAA = unique(GSAA);

totalsizes = GSAA;

GCAA = abs(tGratingDiameter)./dGratingDiameter;
uGCAA = unique(GCAA);

totalcontrasts = GSAA;

untGratingDiameter = unique(tGratingDiameterDeg) ;


%Sizediff 
tGratingContrast = abs(tGratingContrast); 
tlt = find(tLeftTrial == -1);
trt = find(tLeftTrial == 1) ;
sizes = zeros(length(tLeftTrial), 1) ;
sizes(tlt) = dGratingContrast(tlt)./tGratingContrast(tlt);
sizes(trt) = tGratingContrast(trt)./dGratingContrast(trt);
sizediff = sizes;

sizediff = round(sizediff,3);

unsize_plot = unique(sizediff);

%Contrastdiff
tGratingDiameter = abs(tGratingDiameter) ;
tlt = find(tLeftTrial == -1);
trt = find(tLeftTrial == 1) ;
contrasts = zeros(length(tLeftTrial), 1) ;
contrasts(tlt) = dGratingDiameter(tlt)./tGratingDiameter(tlt);
contrasts(trt) = tGratingDiameter(trt)./dGratingDiameter(trt);
contrastdiff = contrasts;



sizediff = round(sizediff, 3) ;
contrastdiff = round(contrastdiff, 3) ;

unsizediff = unique(sizediff);
uncontrastdiff = unique(contrastdiff);

for i = 1:length(uncontrastdiff) ;
conidx = find(contrastdiff == uncontrastdiff(i));
sizecon = sizediff(conidx); %find sizediff data at each contrastdiff
tRRsizecon = tRightResponse(conidx);
unsizecon = unique(sizecon); %find unique sizes within each contrastdiff
for f = 1:length(unsizecon) ;
sizecontrials = find(sizecon == unsizecon(f))'; %find trials at each unique size diff 
trialnum(f,i) = length(sizecontrials); %find number of trials at each unique size diff 
trialresp(f,i) = length(find(tRRsizecon(sizecontrials) == 1)); %find number of right responses for these trials
end
end

find(contrastdiff == uncontrastdiff(1));
contrastdiff(ans) =  uncontrastdiff(7);
find(contrastdiff == uncontrastdiff(2));
contrastdiff(ans) = uncontrastdiff(6);
find(contrastdiff == uncontrastdiff(3)) ;
contrastdiff(ans) = uncontrastdiff(5);
find(contrastdiff == uncontrastdiff(4));
contrastdiff(ans) = uncontrastdiff(4);

uncontrastdiff = unique(contrastdiff);

for i = 1:length(uncontrastdiff) ;
conidx = find(contrastdiff == uncontrastdiff(i));
sizecon = tGratingDiameterDeg(conidx); 
tSuccsizecon = SIx(conidx);
untGratingContrast = unique(sizecon) ;
for f = 1:length(untGratingContrast) ;
sizecontrials = find(sizecon == untGratingContrast(f))';
trialnum(f,i) = length(sizecontrials);
trialresp(f,i) = length(find(tSuccsizecon(sizecontrials) == 1));
end
end



diffsize_con = trialresp./trialnum;

%colormap grey 
figure;
hold on
for i = 1:length(uncontrastdiff);
[PHAT, PCI] = binofit(trialresp(:,i), trialnum(:,i), .05); 
PCI(:,1) = diffsize_con(:,i)- PCI(:,1);
PCI(:,2) = PCI(:,2) - diffsize_con(:,i);
hold on; errorbar(unsizecon, diffsize_con(:,i), PCI(:,1), PCI(:,2)); 
end

figure;
for i = 1:length(uncontrastdiff)
hold on; errorbar(unsizecon, diffsize_con(:,i), PCI(:,1), PCI(:,2), 'LineWidth', 2);
end

figure;
for i = 1:length(uncontrastdiff)
hold on; plot(unsize_plot, diffsize_con(:,i),'LineWidth', 2);
end


co = [0 0 0;
      .1 .1 .1;
      .3 .3 .3;
      .5 .5 .5;
      .7 .7 .7;
     .9 .9 .9;]
set(groot,'defaultAxesColorOrder',co)


co = [0 0 0;
      .1 .1 .1;
      .2 .2 .2;
      .3 .3 .3;
      .4 .4 .4;
      .5 .5 .5;
      .6 .6 .6;
      .7 .7 .7;
      .8 .8 .8;
      .9 .9 .9;]
set(groot,'defaultAxesColorOrder',co)

co = [.1286 .1286 .1286;
      .2571 .2571 .2571;
      .3857 .3857 .3857;
      .5143 .5143 .5143;
      .6429 .6429 .6429;
      .7715 .7715 .7715;
      .9 .9 .9;]
set(groot,'defaultAxesColorOrder',co)

ax = gca ;
ax.XScale = 'log' 
xlim([0 21]) 

hold off
%% model 

sizediff = tGratingDiameterDeg./dGratingDiameterDeg
contrastdiff = tGratingContrast./dGratingContrast
contrastdiff = round(contrastdiff, 3) 
uncontrastdiff = unique(contrastdiff)

GC = GSAA

unGC = unique(GC) 

consAll = unGC
condiscAll = unique(abs(tGratingDiameterDeg))

con1All = zeros(length(GC),1)
oneAll = find(GC(:,1) == consAll(1)); 
con1All(oneAll) = 1
t1All = zeros(length(GC),1) 
toneAll = find(abs(tGratingDiameterDeg(:,1)) == condiscAll(1))
t1All(toneAll) = 1

con2All = zeros(length(GC),1) 
twoAll = find(GC(:,1) == consAll(2)) 
con2All(twoAll) = 1
t2All = zeros(length(GC),1) 
ttwoAll = find(abs(tGratingDiameterDeg(:,1)) == condiscAll(2))
t2All(ttwoAll) = 1 

con3All = zeros(length(GC),1) 
threeAll = find(GC(:,1) == consAll(3)) 
con3All(threeAll) = 1
t3All = zeros(length(GC),1) 
tthreeAll = find(abs(tGratingDiameterDeg(:,1)) == condiscAll(3)) 
t3All(tthreeAll) = 1

con4All = zeros(length(GC),1) 
fourAll = find(GC(:,1) == consAll(4)) 
con4All(fourAll) = 1
t4All = zeros(length(GC),1) 
tfourAll = find(abs(tGratingDiameterDeg(:,1)) == condiscAll(4)) 
t4All(tfourAll) = 1

con1 = zeros(length(GC),1) 
con1(find(contrastdiff == uncontrastdiff(1))) = 1

con2 = zeros(length(GC),1) 
con2(find(contrastdiff == uncontrastdiff(2))) = 1 

con3 = zeros(length(GC),1) 
con3(find(contrastdiff == uncontrastdiff(3))) = 1 

con4 = zeros(length(GC),1) 
con4(find(contrastdiff == uncontrastdiff(4))) = 1

con5 = zeros(length(GC),1) 
con5(find(contrastdiff == uncontrastdiff(5))) = 1

con6 = zeros(length(GC),1) 
con6(find(contrastdiff == uncontrastdiff(6))) = 1

con7 = zeros(length(GC),1) 
con7(find(contrastdiff == uncontrastdiff(7))) = 1

con8 = zeros(length(GC),1) 
con8(find(contrastdiff == uncontrastdiff(8))) = 1

con9 = zeros(length(GC),1) 
con9(find(contrastdiff == uncontrastdiff(9))) = 1



interceptAll = ones(length(tLeftTrial),1)
table_size_contrast = table(SIx, con1All, con2All, con3All, con4All, t1All, t2All, t3All, t4All, con1, con2, con3, con4, con5, con6, con7, con8, con9, tLeftTrial) 
model_size_contrast = 'SIx ~ con1All + con2All + con3All + con4All + t1All + t2All + t3All + t4All + con1 + con2 + con3 + con4 + con5 + con6 + con7 + con8 + con9 + tLeftTrial - 1'

glm_size_contrast = fitglm(table_size_contrast, model_size_contrast, 'distribution', 'binomial') 

glm_size_contrast_estimates = glm_size_contrast.Coefficients.Estimate


figure; hold on 
subplot(2,2,1)
plot(unGC, glm_size_contrast_estimates(1:4),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 40])
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Size Ratio (R/L)')

hold on
subplot(2,2,2)
plot(untGratingDiameter(5:8), glm_size_contrast_estimates(5:8),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 40])
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Target Size')

hold on
subplot(2,2,3)
axis('square')
plot(uncontrastdiff, glm_size_contrast_estimates(9:17),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 5]) 
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Contrast Ratio (T/D)')


hold on
subplot(2,2,4)
axis('square')
plot(1, glm_size_contrast_estimates(18),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 2]) 
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xticks([1])
xticklabels(['Bias'])

untGratingContrast = unique(tGratingContrast)

con1 = zeros(length(GC),1) 
con1(find(tGratingContrast == untGratingContrast(1))) = 1

con2 = zeros(length(GC),1) 
con2(find(tGratingContrast == untGratingContrast(2))) = 1

con3 = zeros(length(GC),1) 
con3(find(tGratingContrast == untGratingContrast(3))) = 1

con4 = zeros(length(GC),1) 
con4(find(tGratingContrast == untGratingContrast(4))) = 1

con5 = zeros(length(GC),1) 
con5(find(tGratingContrast == untGratingContrast(5))) = 1

table_size_contrast_targets = table(SIx, con1All, con2All, con3All, con4All, t1All, t2All, t3All, t4All, con1, con2, con3, con4, con5, tLeftTrial) 
model_size_contrast_targets = 'SIx ~ con1All + con2All + con3All + con4All + t1All + t2All + t3All + t4All + con1 + con2 + con3 + con4 + con5 + tLeftTrial - 1'

glm_size_contrast_targets = fitglm(table_size_contrast_targets, model_size_contrast_targets, 'distribution', 'binomial') 

glm_size_contrast_estimates_targets = glm_size_contrast_targets.Coefficients.Estimate


figure; hold on 
subplot(2,2,1)
plot(unGC, glm_size_contrast_estimates_targets(1:4),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 40])
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Size Ratio (R/L)')

hold on
subplot(2,2,2)
plot(untGratingDiameter(5:8), glm_size_contrast_estimates_targets(5:8),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 40])
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Target Size')

hold on
subplot(2,2,3)
axis('square')
plot(untGratingContrast, glm_size_contrast_estimates_targets(9:13),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 5]) 
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xlabel('Target Contrast')


hold on
subplot(2,2,4)
axis('square')
plot(1, glm_size_contrast_estimates_targets(14),'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
ylim([-5 5])
xlim([0 2]) 
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'
xticks([1])
xticklabels(['Bias'])


%% plot summary 

figure; 
hold on 
xlim([0 1]) 
ylim([-1 6])
axis('square')
plot(i593_target_contrasts, i593_glm_size_contrast_estimates_targets(9:13), 'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
plot(i600_target_contrasts, i600_glm_size_contrast_estimates_targets(9:12), 'o', 'MarkerEdgeColor', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]) 
h = refline(0,0)
h.LineStyle = '--'
h.Color = 'k'

all_targcon_estimates = horzcat(i593_glm_size_contrast_estimates_targets(10:13), i600_glm_size_contrast_estimates_targets(9:12))
for i = 1:length(all_targcon_estimates) 
all_targcon_mean(i) = mean(all_targcon_estimates(i,:))
end
for i = 1:length(all_targcon_estimates) 
all_targcon_sem(i) = std(all_targcon_estimates(i,:))./sqrt(length(all_targcon_estimates(i,:)'))
end
errorbar(i600_target_contrasts, all_targcon_mean, all_targcon_sem, 'ko', 'MarkerFaceColor', 'k') 

cd('G:\home\lindsey\Manuscripts\PerceivedContrast\Figure4')
print('figure4g_targetcontrastweights','-dpdf')
savefig('figure4g_targetcontrastweights')


