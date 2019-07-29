
%Matrices of important metrics on each day
SuccessesPerDay = [];
IncorrectsPerDay = [];
TrialsPerDay = [];


for i=1:length(days_to_include)
    SuccessesPerDay = [SuccessesPerDay sum(s(i).SIx)];
    IncorrectsPerDay = [IncorrectsPerDay sum(s(i).FIx)];
    TrialsPerDay = [TrialsPerDay length(s(i).SIx)];
end

%Combine data from different days into single matrices. 
successes = [];
incorrects = [];
ignores = [];
lefttrials = [];
DecisionTimeMs = []; 
ori = [];


for i = 1:length(days_to_include)
    successes = [successes s(i).SIx];
    incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
    ignores = [ignores s(i).tIgnore];
    lefttrials = [lefttrials s(i).tLeftTrial];
    ori = [ori s(i).tGratingDirectionStart];
    DecisionTimeMs = [DecisionTimeMs s(i).tDecisionTimeMs];
end

DecisionTimeS = DecisionTimeMs./1000;
TrialNum = [1:length(DecisionTimeS)];

figure; scatter(TrialNum, DecisionTimeS)
hold on
xlim([0,2300])
title({'Reaction Time as a Function of Trial Number in Early Training Stages';'Subject: i1401';'Dates: 03/27-04/04'})
xlabel('Trial Number')
ylabel('Reaction Time (s)')
hold off

%%
%Training Day vs y1: Diff Rxn time Success-incorrect and vs y2: Selectivity
 
MeanIncorrect = [];
MeanCorrect = [];
for i=1:length(days_to_include)
    MeanIncorrect= [MeanIncorrect median(s(i).tDecisionTimeMs(logical(s(i).FIx)))];
    MeanCorrect = [MeanCorrect median(s(i).tDecisionTimeMs(logical(s(i).SIx)))];
end

RxnTimeDiff = [];
for i=1:length(days_to_include)
    RxnTimeDiff = [RxnTimeDiff mean(s(i).tDecisionTimeMs(logical(s(i).FIx))) - mean(s(i).tDecisionTimeMs(logical(s(i).SIx)))];
end
RxnTimeDiffS = RxnTimeDiff./1000;

allIncorrect = [];
allCorrect = [];
for i=1:length(days_to_include)
    allIncorrect= [allIncorrect (s(i).tDecisionTimeMs(logical(s(i).FIx)))];
    allCorrect = [allCorrect (s(i).tDecisionTimeMs(logical(s(i).SIx)))];
end

allIncorrectMean = median(allIncorrect,2);
allCorrectMean = median(allCorrect,2);

MeanIncorrect = MeanIncorrect./1000;
MeanCorrect = MeanCorrect./1000;
aaDiff = MeanIncorrect-MeanCorrect;

%Code for selectivity
LeftResponses45 = [];
LeftResponsesNeg45 = [];
Trials45 = [];
TrialsNeg45 = [];


for i=1:length(days_to_include)
    LeftResponses45  = [LeftResponses45  sum(s(i).tLeftResponse(s(i).tGratingDirectionStart==45))]; %Need to change tLeftTrial variable to tLeftResponse
        LeftResponsesNeg45 = [LeftResponsesNeg45 sum(s(i).tLeftResponse(s(i).tGratingDirectionStart==-45))];
        Trials45 = [Trials45 sum(s(i).tGratingDirectionStart==45)];
        TrialsNeg45 = [TrialsNeg45 sum(s(i).tGratingDirectionStart==-45)];
end

selectivity = [];
for i=1:length(days_to_include)
    selectivity = [selectivity (LeftResponses45(i)./Trials45(i)) - (LeftResponsesNeg45(i)./TrialsNeg45(i))];
end

figure;
yyaxis left
plot(days_to_include, MeanCorrect,'b-o');
hold on
plot(days_to_include, MeanIncorrect,'g-o');
ylabel('Average Reaction Time (s)')

yyaxis right 
scatter(days_to_include, selectivity,'filled');
ylabel('Selectivity')
ylim([0.7,1])
xlim([1,8])
title({'Average Reaction Time and Selectivity in Late Training Stages';'Subject: i1401';'Dates: 04/23-05/07'})
xlabel('Day of Training')
legend('Correct Trials', 'Incorrect Trials')
hold off

        
        
        
        
 
   
%  LeftResponse = [];
% 
% for i = 1:length(lefttrials)
%     if lefttrials(i)==1 & successes(i)==1
%         LeftResponse(i)=1;
%     elseif lefttrials(i)==1 & successes(i)==0
%         LeftResponse(i)=0;
%     elseif lefttrials(i)==0 & successes(i)==1
%         LeftResponse(i)=0;
%     elseif lefttrials(i)==0 & successes(i)==0
%         LeftResponse(i)=1;
%     end
% end
% 
% %45 and -45 degree trials only 
% PosTrials = [];
% NegTrials = [];
% 
% for i=1:length(ori)
%     if ori(i)==45
%         PosTrials(i)=1;
%     elseif ori(i)~=45
%         PosTrials(i)=0;
%     end
% end
% 
% for i=1:length(ori)
%     if ori(i)==-45
%         NegTrials(i)=1;
%     elseif ori(i)~=-45
%         NegTrials(i)=0;
%     end
% end
% 
% %Check
% Pos45Trials = ori==45;
% sum(Pos45Trials)
% sum(PosTrials)
% Neg45Trials = ori==-45;
% sum(Neg45Trials)
% sum(NegTrials)
% 
% Pos45Trials = logical(Pos45Trials);
% LeftResponse = logical(LeftResponse);
% Pos45Correct = find(Pos45Trials(LeftResponse));
% Pos45pct = length(Pos45Correct)/sum(Pos45Trials)
% 
% Neg45Trials = logical(Neg45Trials);
% RightResponse = ~LeftResponse; %I think this will include ignores, darn
% Neg45LeftResponse = find(Neg45Trials(RightResponse));
% Neg45pctL = length(Neg45LeftResponse)/sum(Neg45Trials)

        

        

