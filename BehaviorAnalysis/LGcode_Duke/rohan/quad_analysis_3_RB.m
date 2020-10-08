% quad multi day concat attempt

cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')

for x=1:length(input)
    tlt_length(x)=length(input{x}.tLeftTrial);
end

max_tlt=nanmax(tlt_length);
tLeftTrial=nan(length(input), max_tlt);
for x=1:length(input)
    clear temp_tlt
    temp_tlt=double(celleqel2mat_padded((input{x}.tLeftTrial)));
    tLeftTrial(x,1:length(temp_tlt)) = temp_tlt;
end
clear tlt_length; clear x; clear temp_tlt

%% repeat TLT matrix read in for quad values
clear num_trials max_trials qvals
for x=1:length(input)
    num_trials(x)=length(input{x}.tLeftTrial);
end
max_trials=nanmax(num_trials);
qvals=nan(length(input),max_trials,100000);
%qvals=nan(length(input),max_trials);
for x=1:length(input)
    for i=1:num_trials(x)
        clear temp_qval
        temp_qval=double(((input{x}.quadratureValues{i})));
        qvals(x,i,1:length(temp_qval)) = temp_qval;
    end
end

%%
for i=1:length(input)
    tGratingContrast(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tGratingContrast));
    dGratingContrast(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.dGratingContrast));
    SIx(i,1:length(input{i}.tLeftTrial))=double(strcmp(input{i}.trialOutcomeCell, 'success'));
    FIx(i,1:length(input{i}.tLeftTrial))=double(strcmp(input{i}.trialOutcomeCell, 'incorrect'));
    tGratingDiameter(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tGratingDiameterDeg));
    dGratingDiameter(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.dGratingDiameterDeg));
    tBlock2(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tBlock2TrialNumber));
    rightresponse(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tRightResponse));
    leftresponse(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tLeftResponse));
    decisiontime(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.tDecisionTimeMs));
    stimonTime(i,1:length(input{i}.tLeftTrial))=double(celleqel2mat_padded(input{i}.stimTimestampMs));
    
end

  
%%
dGratingDiameter=celleqel2mat_padded(input.dGratingDiameterDeg);
tBlock2=celleqel2mat_padded(input.tBlock2TrialNumber);
rightresponse=celleqel2mat_padded(input.tRightResponse);
leftresponse=celleqel2mat_padded(input.tLeftResponse);
decisiontime=celleqel2mat_padded(input.tDecisionTimeMs);
stimonTimeMs=input.stimOnTimeMs;
decisionthreshold=input.leftDecisionThreshold;
ignores=intersect(find(~SIx), find(~FIx));
ignoredtrials=zeros(1,length(tLeftTrial));
ignoredtrials(ignores)=1;
qfast=celleqel2mat_padded(input.qStartReact);
tRightTrial=zeros(1,length(tLeftTrial));
right_trials=(find(~tLeftTrial));
tRightTrial(right_trials)=1;
%%
for x=1:length(input)
    for i=1:num_trials(x)
        temp_stim_time=stimonTimeMs(x,i);
        for c=1:length(qvals(x,i)
            % dont use loop, use find and then take minimum or start value of vector
            % (working with time so first val will always be min)
            
    
