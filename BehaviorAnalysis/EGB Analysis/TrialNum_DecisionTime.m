%Drag single day data file over

Iix = find(strcmp(input.trialOutcomeCell, 'ignore')); %command strcmp says match "trialOutcomeCell" when it says "ignore"
%Find all ignore trials

%Nix = find(celleqel2mat_padded(input.didNoGo)); 

Tix = setdiff(1:length(input.trialOutcomeCell), Iix); %"setdiff" means give 2 inputs, what is in 1 but not the other? 
%Creating an array of not-ignores

StimStart = input.stimTimestampMs;
NewTrialStart = input.tThisTrialStartTimeMs;

LeftTrials = double(cell2mat(input.tLeftTrial));
TLeftTrials = LeftTrials(Tix); %don't need to use find because just indexing 
LeftTrials = TLeftTrials';

%ADDING ABILITY TO DIVIDE BY SUCCESS AND INCORRECT
TrialOutcome = input.trialOutcomeCell;
TTrialOutcome = TrialOutcome(Tix);
nTrials = length(Tix);
TrialOutcomeMat = nan(nTrials,1);

%Create matrix of trialoutcomes where 1=success, 0=incorrect
for itrial = 1:nTrials
    clear temp 
    temp = TTrialOutcome{itrial}; %calling out a cell
    if strcmp(temp, 'success')  %strcmp = string compare- check if identital'' means string
        TrialOutcomeMat(itrial) = 1;
    elseif strcmp(temp, 'incorrect') 
        TrialOutcomeMat(itrial) = 0;
    end 
end

%maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
% %Max decision time; will be single value 

%Input ITI & adapt time to get proper time window before stimOn

itiTime = double(input.itiTimeMs);
if input.doAdapt==1
    AdaptPeriod = double(input.adaptPeriodMs);
elseif input.doAdapt==0
    AdaptPeriod = 0;
end

%Load matrix of decision times
DecisionTimeMs = double(cell2mat(input.tDecisionTimeMs));
DecisionTimeS = DecisionTimeMs./1000;

%Make matrix counting from 1:max trials for the day, including ignores
TrialNum = [1:length(DecisionTimeMs)];

%Plot scatter of TrialNum vs DecisionTime
figure; scatter(TrialNum, DecisionTimeS)
hold on
title({'Reaction Time as a Function of Trial Number in Early Training Stages';'Subject: i1401';'Date: 03/29/19'})
xlabel('Trial Number')
ylabel('Reaction Time (s)')
hold off

% figure; errorbar(Auniqueori,A90pctL,A90CI(1,:),A90CI(2,:),'r','LineWidth',1.5)
% hold on
% ylim([0,1])
% errorbar(Auniqueori,A0pctL,A0CI(1,:),A0CI(2,:),'g','LineWidth',1.5)
% errorbar(Auniqueori,AnapctL,AnaCI(1,:),AnaCI(2,:),'b','LineWidth',1.5)
% title({'Psychometric Function per Adaptation Condition'; 'Subject: i601'; '12/20 - 02/07, n=8599 trials';'Adapt Period: 4 s'})
% xlabel('Target Orientation (degrees)')
% ylabel('% Left Response')
% legend('90 deg adapter, n=3265','0 deg adapter, n=3414','no adapter, n=1920')
% hold off



%%
%Are there ignores on this day?
ignores=[];
ntrials = length(TrialOutcome);

for itrial=1:ntrials
clear temp
temp = TrialOutcome{itrial};
if strcmp(temp,'ignore')
ignores(itrial)=1;
elseif strcmp(temp,'success')
ignores(itrial)=0;
elseif strcmp(temp,'incorrect')
ignores(itrial)=0;
end
end

IgnoreTrials = find(ignores);

