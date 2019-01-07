%Get matrix of rxn times for all trials across days 

RxnTime = s(1).tDecisionTimeMs;
for i = 2:length(s); 
    RxnTime_add = s(i).tDecisionTimeMs; 
    RxnTime = horzcat(RxnTime, RxnTime_add);
end 

%%
%Get matrix of only Left trials

LTrials=find(s(1).tLeftTrial);
for i=2:length(s);
    LTrials_add=find(s(i).tLeftTrial);
    LTrials = horzcat(LTrials,LTrials_add);
end

%%
%Get matrix of rxn times for Left trials only. How know they are correct
%rxn times and not just the first 1544 in the rxn time matrix?

LTrials_idx=zeros(1,length(RxnTime))
LTrials_idx(LTrials)=1;

RxnTimesL=RxnTime(LTrials)
%%
%get matrix of only Right trials

RTrials=find(~LTrials_idx);
RTrials_idx=zeros(1,length(RTrials));
RTrials_idx(RTrials)=1;

RxnTimesR=RxnTime(RTrials);

%%

figure(1)
hold on
cdfplot(RxnTimesL)
cdfplot(RxnTimesR)
cdfplot(RxnTime)
legend('Left' , 'Right', 'Avg')
%%
