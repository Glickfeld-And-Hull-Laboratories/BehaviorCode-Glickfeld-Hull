
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

tGratingContrast = s_cropped_include.tGratingDiameter;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).tGratingDiameter;
    tGratingContrast = horzcat(tGratingContrast, TGC); 
end 

dGratingContrast = s_cropped_include.dGratingDiameter;
for i = 2:length(s_cropped_include); 
    DGC = s_cropped_include(i).dGratingDiameter;
    dGratingContrast = horzcat(dGratingContrast, DGC); 
end 

tRR = s_cropped_include.tRR;
for i = 2:length(s_cropped_include); 
    tR_resp = s_cropped_include(i).tRR;
    tRR = horzcat(tRR, tR_resp); 
end 

tIgnore = s_cropped_include.tIgnore;
for i = 2:length(s_cropped_include); 
    Ignore = s_cropped_include(i).tIgnore;
    tIgnore = horzcat(tIgnore, Ignore); 
end 

tBlock2 = s_cropped_include.tBlock2; 
for i = 2:length(s_cropped_include); 
    TBT = s_cropped_include(i).tBlock2; 
    tBlock2 = horzcat(tBlock2, TBT); 
end 

nTrials = length(tLeftTrial); 
tIgnore = zeros(1,nTrials); 
tIgnore(intersect(find(SIx == 0 ),find(FIx == 0))) = 1;
ignored_trials = find(~tIgnore) 

SIx = SIx' 
FIx = FIx'
tGratingContrast = tGratingContrast' 
dGratingContrast = dGratingContrast'
tBlock2 = tBlock2' 
tIgnore = tIgnore'
tLeftTrial = tLeftTrial'
tRR = tRR'



SIx = SIx(ignored_trials)
FIx = FIx(ignored_trials)
tGratingContrast = tGratingContrast(ignored_trials)
dGratingContrast = dGratingContrast(ignored_trials)
tBlock2 = tBlock2(ignored_trials) 
tIgnore = tIgnore(ignored_trials)
tLeftTrial = tLeftTrial(ignored_trials)
tRR = tRR(ignored_trials)




tGratingContrast = abs(tGratingContrast) 
tlt = find(tLeftTrial == 1)
trt = find(tLeftTrial == 0) 
contrasts = zeros(length(tLeftTrial), 1) 
contrasts(tlt) = dGratingContrast(tlt)./tGratingContrast(tlt)
contrasts(trt) = tGratingContrast(trt)./dGratingContrast(trt)
condiff = contrasts
condiff = round(condiff, 3) 

notNogo = find(condiff ~= 1) 

SIx = SIx(notNogo)
FIx = FIx(notNogo)
tGratingContrast = tGratingContrast(notNogo)
dGratingContrast = dGratingContrast(notNogo)
tBlock2 = tBlock2(notNogo)
tIgnore = tIgnore(notNogo)
tLeftTrial = tLeftTrial(notNogo)
tRR = tRR(notNogo)
condiff = condiff(notNogo) 

x = unique(condiff);

condiff_control = condiff(find(tBlock2 ==0)) 
condiff_stim = condiff(find(tBlock2)) 

tRR_control = tRR(find(tBlock2 == 0))
tRR_stim = tRR(find(tBlock2))


for i = 1:length(x)
numtrials_control(i) = sum(condiff_control == x(i)); 
end
for i = 1:length(x)
numreal_control(i) = length(intersect(find(condiff_control == x(i)), find(tRR_control == 1)));
end

for i = 1:length(x)
numtrials_stim(i) = sum(condiff_stim == x(i)); 
end
for i = 1:length(x)
numreal_stim(i) = length(intersect(find(condiff_stim == x(i)), find(tRR_stim == 1)))
end


numtrials_control = numtrials_control'
numreal_control = numreal_control'


numtrials_stim = numtrials_stim'
numreal_stim = numreal_stim'


[PHATRC_control, PCIRC_control] = binofit(numreal_control, numtrials_control, .05)

[PHATRC_stim, PCIRC_stim] = binofit(numreal_stim, numtrials_stim, .05)

realdatachoice_control(:,1) = numreal_control./numtrials_control

realdatachoice_stim(:,1) = numreal_stim./numtrials_stim


PCIRC_control(:,1) = realdatachoice_control(:,1) - PCIRC_control(:,1)
PCIRC_control(:,2) = PCIRC_control(:,2) - realdatachoice_control(:,1)
PCIRC_stim(:,1) = realdatachoice_stim(:,1) - PCIRC_stim(:,1)
PCIRC_stim(:,2) = PCIRC_stim(:,2) - realdatachoice_stim(:,1)


figure; hold on; 
errorbar(x, realdatachoice_control, PCIRC_control(:,1), PCIRC_control(:,2), 'k') 
errorbar(x, realdatachoice_stim, PCIRC_stim(:,1), PCIRC_stim(:,2), 'r') 
ax = gca 
ax.XScale = 'log' 
xlim([0 101]) 
axis('square')
ylabel('% Right Responses') 
xlabel('Size Ratio (R/L)') 
title('i600 right responses by size ratio') 


