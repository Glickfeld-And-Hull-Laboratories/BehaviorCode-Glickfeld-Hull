%% create a structure for all data
cd('Z:\2018\Code\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis')
for i = 1:length(input)
    s(i).date = input{i,1}.saveTime;
end

for i = 1:length(input)
s(i).SIx = double(strcmp(input{i,1}.trialOutcomeCell,'success'));
end

for i = 1:length(input)
s(i).FIx = double(strcmp(input{i,1}.trialOutcomeCell,'incorrect'));
end

for i = 1:length(input)
s(i).tLeftTrial = celleqel2mat_padded(input{i,1}.tLeftTrial);
end

for i = 1:length(input)
s(i).tGratingDirectionStart = celleqel2mat_padded(input{i,1}.tGratingDirectionStart);
end


for i = 1:length(input)
s(i).tGratingDiameter = celleqel2mat_padded(input{i,1}.tGratingDiameterDeg);
end

for i = 1:length(input)
s(i).tBlock2 = celleqel2mat_padded(input{i,1}.tBlock2TrialNumber);
end

for i = 1:length(input)
s(i).isNoGo = celleqel2mat_padded(input{i,1}.isNoGo);
end

for i = 1:length(input)
s(i).tRR = celleqel2mat_padded(input{i,1}.tRightResponse);
end

for i = 1:length(input)
s(i).tLR = celleqel2mat_padded(input{i,1}.tLeftResponse);
end

for i = 1:length(input)
s(i).ProbLeft = celleqel2mat_padded(input{i,1}.tStimProbAvgLeft);
end

for i = 1:length(s)
s(i).unProbLeft = unique(s(i).ProbLeft);
end


for i = 1:length(s)
s(i).tDecisionTimeMs = celleqel2mat_padded(input{i,1}.tDecisionTimeMs);
end

for i = 1:length(s)
s(i).stimOnTimeMs = input{i,1}.stimOnTimeMs;
end

for i = 1:length(s)
s(i).doFeedbackMotion = input{i,1}.doFeedbackMotion;
end

for i = 1:length(s)
s(i).laserpower = input{i,1}.block2TrialLaserPowerMw;
end

for i = 1:length(s) 
nTrials = length(s(i).tLeftTrial); 
tIgnore = zeros(1,nTrials); 
tIgnore(intersect(find(s(i).SIx == 0 ),find(s(i).FIx == 0))) = 1;
s(i).tIgnore = tIgnore;
end
%% Loop through to plot each day's data, plot ignores, and crop
% define c as trials to crop from earliest date to latest date eg. c = [123 465 222 ...] 
%if there is a day you don't want to include, make c = 0
c = [198, 391, 402, 441, 481, 480, 421, 172, 485, 0, 450, 421, 466, 458, 501, 399, 452, 429, 347, 380, 326, 335]; 

for i = 1:length(s) 
    s(i).SIx = s(i).SIx(1:c(i))
    s(i).FIx = s(i).FIx(1:c(i))
    s(i).tGratingDirectionStart = s(i).tGratingDirectionStart(1:c(i));
    s(i).tGratingDiameter = s(i).tGratingDiameter(1:c(i));
    s(i).tBlock2 = s(i).tBlock2(1:c(i));
    s(i).tLeftTrial = s(i).tLeftTrial(1:c(i));
    s(i).tRR = s(i).tRR(1:c(i));
    s(i).tIgnore = s(i).tIgnore(1:c(i)); 
    s(i).tDecisionTimeMs = s(i).tDecisionTimeMs(1:c(i)); 
end

% if c = 0, s will be an empty cell. Index the full cells and exlude empty.
days_to_include = find(c);
s = s(days_to_include);

%% Exclude days based on selectivity and bias criteria 

for i = 1:length(s)

tLeftTrial = s(i).tLeftTrial;

tGratingDirectionStart = s(i).tGratingDirectionStart;
tGratingDiameter = s(i).tGratingDiameter %was repeated with dGratingDiameter etc for distractor

tRR = s(i).tRR;

tLR = s(i).tLR;

tBlock2 = s(i).tBlock2;

NoBlock = find(tBlock2 == 0);

tLeftTrial = tLeftTrial(NoBlock);
tGratingDirectionStart = tGratingDirectionStart(NoBlock); %was repeated for d
tRR = tRR(NoBlock);
tLR = tLR(NoBlock);
tBlock2 = tBlock2(NoBlock);

tGratingDiameter=tGratingDiameter(NoBlock);%was repeated for d

GC = tGratingDirectionStart./dGratingContrast; 
%GC=tGratingDiameter./dGratingDiameter;
GC = round(GC, 4); 

unGC = unique(GC);


z = length(unGC);

AllResponse = zeros(length(tRR),1); 
AllResponse(find(tRR == 1)) = 1;
AllResponse(find(tLR == 1)) = 1;

GC = GC(find(AllResponse == 1));
tLeftTrial = tLeftTrial(find(AllResponse == 1));
tRR = tRR(find(AllResponse == 1));
tLR = tLR(find(AllResponse == 1));


RightTrials = (intersect(find(GC == unGC(z)), find(tLeftTrial == 0)));
LeftTrials = (intersect(find(GC == unGC(z)), find(tLeftTrial == 1)));

LengthRTrials = length(RightTrials); 
LengthLTrials = length(LeftTrials); 

RightResponse = length(find(tRR(RightTrials) == 1));
LeftResponse = length(find(tLR(LeftTrials) == 1)); 


PercentRight = RightResponse/LengthRTrials;
PercentLeft = 1 - (LeftResponse/LengthLTrials);

Selectivity_cropped(i)  = PercentRight - PercentLeft;

end

cropped_include = find(Selectivity_cropped >= 0); % normally should be above 0.9

s_cropped = s(cropped_include);

input_cropped = input(cropped_include);


for i = 1:length(s_cropped)

tLeftTrial = s_cropped(i).tLeftTrial;

tGratingDirectionStart = s_cropped(i).tGratingContrast;
tGratingDiameter=s_cropped(i).tGratingDiameter;
dGratingContrast = s_cropped(i).dGratingContrast;
dGratingDiameter=s_cropped(i).dGratingDiameter;
tRR = s_cropped(i).tRR;

tLR = s_cropped(i).tLR;

tBlock2 = s_cropped(i).tBlock2;
    
NoBlock = find(tBlock2 == 0);

tLeftTrial = tLeftTrial(NoBlock);
tGratingDirectionStart = tGratingDirectionStart(NoBlock); 
dGratingContrast = dGratingContrast(NoBlock); 
tRR = tRR(NoBlock); 
tLR = tLR(NoBlock); 
tBlock2 = tBlock2(NoBlock); 

tGratingDiameter=tGratingDiameter(NoBlock);
dGratingDiameter=dGratingDiameter(NoBlock); 

GC = tGratingDirectionStart./dGratingContrast;
%GC= tGratingDiameter./dGratingDiameter;
GC = round(GC, 4);

unGC = unique(GC);

AllResponse = zeros(length(tRR),1); 
AllResponse(find(tRR == 1)) = 1;
AllResponse(find(tLR == 1)) = 1;

GC = GC(find(AllResponse == 1)); 
tLeftTrial = tLeftTrial(find(AllResponse == 1));
tRR = tRR(find(AllResponse == 1));
tLR = tLR(find(AllResponse == 1));

for z = 1:length(unGC);
RightTrials = (intersect(find(GC == unGC(z)), find(tLeftTrial == 0)));
LeftTrials = (intersect(find(GC == unGC(z)), find(tLeftTrial == 1)));

LengthRTrials(z) = length(RightTrials);
LengthLTrials(z) = length(LeftTrials); 

RightResponse(z) = length(find(tRR(RightTrials) == 1));
LeftResponse(z) = length(find(tLR(LeftTrials) == 1)); 


PercentRight(z) = RightResponse(z)./LengthRTrials(z); 
PercentLeft(z) = 1 - (LeftResponse(z)./LengthLTrials(z));

N_trials_con(z) = LengthRTrials(z) + LengthLTrials(z);

end

for b = 1:length(unGC)
Bias(b)=((PercentRight(b)+PercentLeft(b))/2-0.5).*(N_trials_con(b)/sum(N_trials_con));
end

B = nansum(Bias);

totalbias1(i) = B;
end

Bmeasure = find(abs(totalbias1) <= 1); %should normally be 0.1

s_cropped_include = s_cropped(Bmeasure);
input_cropped_include = input_cropped(Bmeasure);
 
