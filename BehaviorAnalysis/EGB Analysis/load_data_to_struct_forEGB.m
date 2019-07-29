%% create a structure for all data
%cd('\\crash.dhe.duke.edu\data\home\emily\2018\Code\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis')
cd('C:\Users\emily2\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis');

for i = 1:length(input)
    s(i).subjectNum = input{i,1}.subjectNum;
end

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
s(i).tLeftResponse = celleqel2mat_padded(input{i,1}.tLeftResponse);
end

for i = 1:length(input)
s(i).tGratingDirectionStart = celleqel2mat_padded(input{i,1}.tGratingDirectionStart);
end

% ONLY COMMENT OUT FOR EARLY 601/602
for i = 1:length(input)
s(i).doAdapt = input{i,1}.doAdapt;
end

% ONLY COMMENT OUT FOR EARLY 601/602
for i = 1:length(input)
s(i).aGratingDirectionDeg = celleqel2mat_padded(input{i,1}.aGratingDirectionDeg);
end

for i = 1:length(input)
s(i).tGratingDirectionStart = celleqel2mat_padded(input{i,1}.tGratingDirectionStart);
end

% ONLY COMMENT OUT FOR EARLY 601/602
for i = 1:length(input)
s(i).aGratingContrast = celleqel2mat_padded(input{i,1}.aGratingContrast);
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

%%
%Input lower & upper bounds for each day's data

%Automate this? Crop beginning and end of trials within each day - do this first and then
%day selection?? 

%4sec adapt i601, i1400
% c= [400 531 329 313 498 372 0 448 336 0 354 350 577 0 0 322 291 0 442 438 ...
%     0 440 340 361 366 250 0 322 331 370 303 358 362 0 329 400 363 0 414 369 ...
%     308 488 200 385 300 250 ... %end601 start1400
%     0 229 225 411 513 0 0 150 343 369 296 355 309 339 0 426 350 0 338 381 ...
%     300 200 250 334 250 200 438 200 300 150 250 0];
% b= [1 1 1 1 1 1 1 1 1 1 1 1 500 1 1 1 1 1 1 250 ...
%     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
%     1 100 1 200 1 1 ... %end601 start1400
%     1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
%     1 1 1 1 1 1 1 1 1 1 1 1];

%100ms Adapter 601,1400,1401 
% c = [489 150 300 250 0 0 425 458 481 430 100 200 ... %end601 start1400
%     446 250 400 300 400 412 433 200 421 428 ... %end1400 start1401
%     300 447 422 505 364 445 438 300];
% b = [1 1 1 1 1 1 1 1  1 1 1 1 ... %end601 start 1400
%     1 1 1 1 1 1 1 1 1 300 ... %end1400 start1401
%     1 1 1 1 1 1 1 1];

%Flash Adapt 1401 
% c= [359 409 150 300 300 448 339 416 200 381 389 383 0 0 381 300 435 250 100 0 295 0 375 0 416 419 200];
% b= [1 1 1 1 1 1 1 1 1 1 1 1 1 1 281 200 335 1 1 1 150 1 275 0 150 1 1];
% b= ones(1,length(c));

%Flash Adapt 1401 + i1402 
c = [359 409 150 300 300 448 339 416 200 381 389 383 0 0 381 300 435 250 100 0 295 0 375 0 416 419 200 ...
    200 384 150 0 175 409 390 466 476 400 465 400 384];
b = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 281 200 335 1 1 1 150 1 275 0 150 1 1 ...
    1 1 1 1 1 1 1 1 1 1 1 200 1];

for i = 1:length(s) 
    s(i).SIx = s(i).SIx(b(i):c(i));
    s(i).FIx = s(i).FIx(b(i):c(i));
    s(i).tGratingDirectionStart = s(i).tGratingDirectionStart(b(i):c(i));
    s(i).aGratingDirectionDeg = s(i).aGratingDirectionDeg(b(i):c(i));
    s(i).aGratingContrast = s(i).aGratingContrast(b(i):c(i));
    s(i).tBlock2 = s(i).tBlock2(b(i):c(i));
    s(i).tLeftTrial = s(i).tLeftTrial(b(i):c(i));
    s(i).tLeftResponse = s(i).tLeftResponse(b(i):c(i));
    s(i).tRR = s(i).tRR(b(i):c(i));
    s(i).tIgnore = s(i).tIgnore(b(i):c(i)); 
    s(i).tDecisionTimeMs = s(i).tDecisionTimeMs(b(i):c(i)); 
end

%%
% Exclude days based on selectivity
for i=1:length(s)
    tLeftTrial = s(i).tLeftTrial;
    Success = s(i).SIx(s(i).aGratingContrast==0);
    Failure = s(i).FIx(s(i).aGratingContrast==0);
    N_EasyTrials(i) = sum(abs(s(i).tGratingDirectionStart(s(i).aGratingContrast==0))==45); 
    N_SIx(i) =  sum(Success(abs(s(i).tGratingDirectionStart(s(i).aGratingContrast==0))==45)); 
    N_FIx(i) =  sum(Failure(abs(s(i).tGratingDirectionStart(s(i).aGratingContrast==0))==45)); 
    Selec(i) =  N_SIx(i) ./( N_SIx(i) +N_FIx(i)); 
    if Selec(i)>= 0.90
        c(i)=c(i);
    elseif Selec(i)< 0.90
        c(i)=0;
    end
end

%%
%Exclude days based on bias criteria 
%Single output of bias measure per each i 

for i=1:length(s)
    LRDiff = [];
    OriWeights = [];
    UniqueOri = unique(abs(s(i).tGratingDirectionStart));
    tLeftTrial = s(i).tLeftTrial;
    for ori=1:length(UniqueOri)
        match_grating = abs(s(i).tGratingDirectionStart)==UniqueOri(ori);
        LCorrect = sum(s(i).tLeftTrial .* s(i).SIx .* match_grating);
        LTotal = sum(s(i).tLeftTrial .* s(i).SIx .* match_grating) + sum(s(i).tLeftTrial .* s(i).FIx .* match_grating);
        LPct = LCorrect./LTotal;
        RTrials = ~s(i).tLeftTrial;
        RCorrect = sum(RTrials .* s(i).SIx .* match_grating);
        RTotal = sum(RTrials .* s(i).SIx .* match_grating) + sum(RTrials .* s(i).FIx .* match_grating);
        RPct = RCorrect./RTotal;
        LRDiff = [LRDiff LPct-RPct];
        Nori = sum(s(i).tGratingDirectionStart==UniqueOri(ori));
        NTotal = length(s(i).tGratingDirectionStart);
        OriWeights = [OriWeights Nori/NTotal];
        Bias = abs(sum(LRDiff.*OriWeights));
    end
    if Bias<=0.1
        c(i)=c(i);
    elseif Bias> 0.1
        c(i)=0;
    end
end

e=ones(1,(length(b)));
for i=1:length(c)
    if c(i)==0
        d(i)=0;
    elseif c(i)>0
        d(i) = c(i)-b(i)+1;
    end
end


%% Loop through to plot each day's data, plot ignores, and crop
% define c as trials to crop from earliest date to latest date eg. c = [123 465 222 ...] 
%if there is a day you don't want to include, make c = 0

for i = 1:length(s) 
    s(i).SIx = s(i).SIx(e(i):d(i));
    s(i).FIx = s(i).FIx(e(i):d(i));
    s(i).tGratingDirectionStart = s(i).tGratingDirectionStart(e(i):d(i));
    s(i).aGratingDirectionDeg = s(i).aGratingDirectionDeg(e(i):d(i));
    s(i).aGratingContrast = s(i).aGratingContrast(e(i):d(i));
    s(i).tBlock2 = s(i).tBlock2(e(i):d(i));
    s(i).tLeftTrial = s(i).tLeftTrial(e(i):d(i));
    s(i).tLeftResponse = s(i).tLeftResponse(e(i):d(i));
    s(i).tRR = s(i).tRR(e(i):d(i));
    s(i).tIgnore = s(i).tIgnore(e(i):d(i)); 
    s(i).tDecisionTimeMs = s(i).tDecisionTimeMs(e(i):d(i)); 
end

% for i = 1:length(s) 
%     s(i).SIx = s(i).SIx(1:c(i));
%     s(i).FIx = s(i).FIx(1:c(i));
%     s(i).tGratingDirectionStart = s(i).tGratingDirectionStart(1:c(i));
%     s(i).aGratingDirectionDeg = s(i).aGratingDirectionDeg(1:c(i));
%     s(i).aGratingContrast = s(i).aGratingContrast(1:c(i));
%     s(i).tBlock2 = s(i).tBlock2(1:c(i));
%     s(i).tLeftTrial = s(i).tLeftTrial(1:c(i));
%     s(i).tLeftResponse = s(i).tLeftResponse(1:c(i));
%     s(i).tRR = s(i).tRR(1:c(i));
%     s(i).tIgnore = s(i).tIgnore(1:c(i)); 
%     s(i).tDecisionTimeMs = s(i).tDecisionTimeMs(1:c(i)); 
% end

% if c = 0, s will be an empty cell. Index the full cells and exlude empty.
%Need to do this step differently to include lower limit (b) ?? 
days_to_include = find(c);
s = s(days_to_include);
clearvars -except s input days_to_include;

%1401 1st 4 days of Flash Stim Adapt
%c = [359 409 150 300];

%601 1st 7 days 09/20-10/01/18: c = [189 222 264 225 268 196 209 233];
%1400 1st 7 days 11/21-11/30/18: c = [270 267 468 357 392 283 248];
%602 1st 7 days 9/20 - 9/28/18: c = [198 391 402 441 481 480 421];
%1401 03/27-04/04(479) - 04/11/19: c= [351 309 527 453 519 380 479];
%1402 late training, 04/23-05/07: c = [453 0 200 447 422 505 364 334 438];
%1401 03/27-04/11/19: c= [351 309 527 453 519 380 479 503 460 408 484 541];

%1401 flash stim 05/09-05/23
% b= ones(1,length(c));
% c= [359 409 150 300 300 448 339 416 200 381 389];

%other data sets: 
%601 4 sec adapt: c = [322 291 525 442 438 0 440 340 361 366 250 0 322 331 370 303 358 362 0 0 400 363 0 350 369 308 488 0 0 300 200];

%601 100ms Adapt: c = [489 150 300 250 0 0 425 400 481 430 100 200];

%1400 100ms Adapt 02/08-02/21: c = [446 250 400 300 400 412 433 200 421 428];

%1401 100ms Adapt 04/23-5/07: 
% b= [1 1 1 1 1 1 1 1 1];

%1400 12/27 - 2/20
%c = [513 218 318 264 343 369 296 355 309 300 381 426 400 0 338 381 0 200 275 334 300 200 438 200 300 200 250 0 446 250 400 300 400 412 433 200 421];
%1400 12/18-01/10: c = [308 348 229 225 411 513 218 318 264 343 369 296 355 309 300];

%601, 1400 & 1401 100ms Adapt Combined
% b = ones(1,31);
% c=[489 150 300 250 0 0 425 400 481 430 100 200 446 250 400 300 400 412 433 200 421 428 453 0 300 447 422 505 364 445 438];

