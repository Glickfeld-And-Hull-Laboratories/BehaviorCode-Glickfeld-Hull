%get SIx, get tLeftResponse and tRightResponse
% find quad value at decision time?
% compare decision direction to toofast time quad threshold value
cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')


SIx=double(strcmp(input.trialOutcomeCell, 'success'));
FIx=double(strcmp(input.trialOutcomeCell, 'incorrect'));
tLeftTrial=celleqel2mat_padded(input.tLeftTrial);
tGratingContrast=celleqel2mat_padded(input.tGratingContrast);
dGratingContrast=celleqel2mat_padded(input.dGratingContrast);
tGratingDiameter=celleqel2mat_padded(input.tGratingDiameterDeg);
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
clear right_trials; clear ignores
%% quad value inputs

for x=1:length(tLeftTrial)
     q_length(x)=length(input.quadratureTimesUs{x});
end

max_q=nanmax(q_length);

for x=1:length(tLeftTrial)
    q_length(x)=length(input.quadratureValues{x});
end
max_qval=nanmax(q_length);
qvals=nan(length(tLeftTrial), max_q);
for x=1:length(tLeftTrial)
    clear temp_qval
    temp_qval=double((input.quadratureValues{x}));
    qvals(x,1:length(temp_qval)) = temp_qval;
end
clear temp_qval; clear x 

% quad time inputs
for x=1:length(tLeftTrial)
    q_length(x)=length(input.quadratureTimesUs{x});
end
max_qval=nanmax(q_length);
qtimes=nan(length(tLeftTrial), max_q);
for x=1:length(tLeftTrial)
    clear temp_qval
    temp_qtime=double((input.quadratureTimesUs{x}));
    qtimes(x,1:length(temp_qtime)) = temp_qtime;
end
qtimes=qtimes/1000;
clear temp_qtime; clear x 

%%
correct_trials=find(SIx);
correct_trials_right_idx=zeros(1,length(tLeftTrial));
correct_trials_right=intersect(find(SIx), find(tRightTrial))
correct_trials_right_idx(correct_trials_right)=1;
correct_trials_left_idx=zeros(1,length(tLeftTrial));
correct_trials_left=intersect(find(SIx), find(tLeftTrial))
correct_trials_left_idx(correct_trials_left)=1;
%%
responsevector=zeros(1,length(tLeftTrial)); % value of 1 in response vector indicates matching qfast and decision
for i=1:length(correct_trials_right)
    c=correct_trials_right(i);
    if qfast(c)<(decisionthreshold*1.5)
        responsevector(c)=-1;
    else 
        responsevector(c)=0;
    end
end

for i=1:length(correct_trials_left)
    c=correct_trials_left(i);
    if qfast(c)>(decisionthreshold*1.5)
        responsevector(c)=1;
    else 
        responsevector(c)=0;
    end
end

num_match_right=length(find(responsevector==-1))
num_match_left=length(find(responsevector==1))
percent_match_right=num_match_right/sum(tRightTrial)
percent_match_left=num_match_left/sum(tLeftTrial)
x=1:2
figure;
%% unrelated
figure(1);
hold on
p=plot(qVals_final(:,1:10));
for i=1:10

    set(p(i), 'color', [0.1*i 1/i i*0.05 ]);
end





        
        
    

