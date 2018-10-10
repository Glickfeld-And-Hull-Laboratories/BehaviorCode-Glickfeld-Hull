function [SIx, FIx, tLeftTrial, tGratingContrast, dGratingContrast, tRightResponse, tIgnore, tBlock2, GCAA, unGCAA, target_contrasts, totalcontrasts] = concat_model_data(s_cropped_include) 

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
    DGC = s_cropped_include(i).dGratingContrast;
    dGratingContrast = horzcat(dGratingContrast, DGC); 
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

tBlock2 = s_cropped_include.tBlock2; 
for i = 2:length(s_cropped_include); 
    TBT = s_cropped_include(i).tBlock2; 
    tBlock2 = horzcat(tBlock2, TBT); 
end 

nTrials = length(tLeftTrial); 
tIgnore = zeros(1,nTrials); 
tIgnore(intersect(find(SIx == 0 ),find(FIx == 0))) = 1;
ignored_trials = find(~tIgnore); 

tGratingContrast(find(tLeftTrial)) = tGratingContrast(find(tLeftTrial)).*-1;

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

left_trials = find(tLeftTrial == 1); 
tLeftTrial(left_trials) = -1 ;
right_trials = find(~tLeftTrial) ;
tLeftTrial(right_trials) = 1 ;


SIx = SIx' ;
FIx = FIx';
tGratingContrast = tGratingContrast';
SIx_prev = SIx_prev';
FIx_prev = FIx_prev' ;
dGratingContrast = dGratingContrast';
tBlock2 = tBlock2' ;
tIgnore = tIgnore';
tLeftTrial = tLeftTrial';
tRightResponse = tRightResponse';


SIx = SIx(ignored_trials);
FIx = FIx(ignored_trials);
tGratingContrast = tGratingContrast(ignored_trials);
SIx_prev = SIx_prev(ignored_trials);
FIx_prev = FIx_prev(ignored_trials);
dGratingContrast = dGratingContrast(ignored_trials);
tBlock2 = tBlock2(ignored_trials);
tLeftTrial = tLeftTrial(ignored_trials);
tRightResponse = tRightResponse(ignored_trials);

GCAA = abs(tGratingContrast)./dGratingContrast;
GCAA = round(GCAA, 3);
uGCAA = unique(GCAA);

GCNG = GCAA;
notNogo = find(GCNG ~= 1);

SIx = SIx(notNogo);
FIx = FIx(notNogo); 
tGratingContrast = tGratingContrast(notNogo); 
tLeftTrial = tLeftTrial(notNogo);
SIx_prev = SIx_prev(notNogo);
FIx_prev = FIx_prev(notNogo);
tIgnore = tIgnore(notNogo); 
dGratingContrast = dGratingContrast(notNogo);
tBlock2 = tBlock2(notNogo);
GCAA = GCAA(notNogo);
unGCAA = unique(GCAA);
tRightResponse = tRightResponse(notNogo)

target_contrasts = unique(tGratingContrast); 

tGratingContrast = abs(tGratingContrast);
tlt = find(tLeftTrial == -1);
trt = find(tLeftTrial == 1); 
contrasts = zeros(length(tLeftTrial), 1); 
contrasts(tlt) = dGratingContrast(tlt)./tGratingContrast(tlt);
contrasts(trt) = tGratingContrast(trt)./dGratingContrast(trt);
contrasts = round(contrasts, 3); 

totalcontrasts = contrasts;

if length(unique(totalcontrasts)) > 8 && length(unique(tGratingContrast)) == 4;
    

d1 = find(totalcontrasts <= .01 & totalcontrasts > 0);
d2 = find(totalcontrasts <= .3 & totalcontrasts > .01);
d3 = find(totalcontrasts <= .5 & totalcontrasts > .3); 
d4 = find(totalcontrasts <= 1 & totalcontrasts > .5); 
d5 = find(totalcontrasts <= 2 & totalcontrasts > 1); 
d6 = find(totalcontrasts <= 5 & totalcontrasts > 2); 
d7 = find(totalcontrasts <= 16 & totalcontrasts > 5); 
d8 = find(totalcontrasts <= 102 & totalcontrasts > 16); 

GCAAU = unique(totalcontrasts); 

du1 = find(GCAAU <= .01 & GCAAU > 0);
du2 = find(GCAAU <= .3 & GCAAU > .01);
du3 = find(GCAAU <= .5 & GCAAU > .3); 
du4 = find(GCAAU <= 1 & GCAAU > .5); 
du5 = find(GCAAU <= 2 & GCAAU > 1); 
du6 = find(GCAAU <= 5 & GCAAU > 2); 
du7 = find(GCAAU <= 16 & GCAAU > 5); 
du8 = find(GCAAU <= 102 & GCAAU > 16); 

totalcontrasts(d1) = GCAAU(du1(end));
totalcontrasts(d2) = GCAAU(du2(end));
totalcontrasts(d3) = GCAAU(du3(end));
totalcontrasts(d4) = GCAAU(du4(end));
totalcontrasts(d5) = GCAAU(du5(1));
totalcontrasts(d6) = GCAAU(du6(1));
totalcontrasts(d7) = GCAAU(du7(1));
totalcontrasts(d8) = GCAAU(du8(1));

elseif length(unique(totalcontrasts)) > 8 && length(unique(tGratingContrast)) == 3;
d1 = find(totalcontrasts <= .01 & totalcontrasts > 0);
d2 = find(totalcontrasts <= .5 & totalcontrasts > .01); 
d3 = find(totalcontrasts <= 1 & totalcontrasts > .5); 
d4 = find(totalcontrasts <= 2 & totalcontrasts > 1); 
d5 = find(totalcontrasts <= 16 & totalcontrasts > 2); 
d6 = find(totalcontrasts <= 102 & totalcontrasts > 16); 

GCAAU = unique(totalcontrasts); 

du1 = find(GCAAU <= .01 & GCAAU > 0);
du2 = find(GCAAU <= .5 & GCAAU > .01); 
du3 = find(GCAAU <= 1 & GCAAU > .5); 
du4 = find(GCAAU <= 2 & GCAAU > 1); 
du5 = find(GCAAU <= 16 & GCAAU > 2); 
du6 = find(GCAAU <= 102 & GCAAU > 16); 

totalcontrasts(d1) = GCAAU(du1(end));
totalcontrasts(d2) = GCAAU(du2(end));
totalcontrasts(d3) = GCAAU(du3(end));
totalcontrasts(d4) = GCAAU(du4(1));
totalcontrasts(d5) = GCAAU(du5(1));
totalcontrasts(d6) = GCAAU(du6(1));

end 
end
