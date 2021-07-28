%% make 2 part subplot of quadrature correct and incorrect trials 

cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
savedate=input.saveTime;
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
stimtime_vector=input.stimTimestampMs;
SIx = double(strcmp(input.trialOutcomeCell,'success'));
FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
block2=celleqel2mat_padded(input.tBlock2TrialNumber);
tContrast=celleqel2mat_padded(input.tGratingContrast); % comment out for size discrim
dContrast=celleqel2mat_padded(input.dGratingContrast); % comment out for size discrim
%tContrast=celleqel2mat_padded(input.tGratingDiameterDeg);
%dContrast=celleqel2mat_padded(input.dGratingDiameterDeg);
contrastratio=tContrast./dContrast;
contrastratio=round(contrastratio,3);
uncontrasts=unique(contrastratio);
c1=uncontrasts(4); c2=uncontrasts(3); c3=uncontrasts(2); c4=uncontrasts(1);
tLeftTrial=celleqel2mat_padded(input.tLeftTrial);
numtrials=length(tLeftTrial);
stimtime_vector=celleqel2mat_padded(input.stimTimestampMs);
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        stimTime = double(input.stimTimestampMs{trN});
        qTimes_zero = qTimes-stimTime;
        qVals = qVals-qVals(1);
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000);
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind);
            qVals_sub = qVals(time_ind);
            qTimes_temp = qTimes(time_ind);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000;
            qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)';
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)';
            if input.tDecisionTimeMs{trN} < 10000
                if isnan(qVals_final(8000,trN))
                    qVals_final(8000,trN) = qVals_final(find(~isnan(qVals_final(:,trN)),1,'first'),trN);
                end
                qVal_off = qVals_final(:,trN)-qVals_final(8000,trN);
                qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN);
            end
        else
            return
        end
    end
end
%%
nledcond=2;
ncorrectcond=2;
nstimlocation=2;
nratios=4;
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1
sorted_trials = nan(nledcond, ncorrectcond, nstimlocation, nratios, numtrials);
for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear tempcon tempbehavesuccess tempbehavefailure tempbehave temploc temprat temp
               tempcon = find(block2 == icontrol-1);
               if ibehavior == 1
                   tempbehave = find(FIx == 1);
               elseif ibehavior == 2
                   tempbehave = find(SIx == 1);
               end
               temploc = find(tLeftTrial == iloc-1);
               temprat = find(contrastratio == uncontrasts(irat));
               temp =intersect(intersect(intersect(tempcon, tempbehave), temploc), temprat);
               sorted_trials(icontrol, ibehavior, iloc, irat, 1:length(temp)) = temp;
           end
       end
   end
end
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1
control_right_correctc1 = squeeze(sorted_trials(1,2,1,4,:));
control_right_correctc2= squeeze(sorted_trials(1,2,1,3,:));
control_right_correctc3= squeeze(sorted_trials(1,2,1,2,:));
control_right_correctc4= squeeze(sorted_trials(1,2,1,1,:));
control_left_correctc1= squeeze(sorted_trials(1,2,2,4,:));
control_left_correctc2= squeeze(sorted_trials(1,2,2,3,:));                
control_left_correctc3= squeeze(sorted_trials(1,2,2,2,:));  
control_left_correctc4= squeeze(sorted_trials(1,2,2,1,:));
control_right_incorrectc1= squeeze(sorted_trials(1,1,1,4));
control_right_incorrectc2= squeeze(sorted_trials(1,1,1,3));
control_right_incorrectc3= squeeze(sorted_trials(1,1,1,2));
control_right_incorrectc4= squeeze(sorted_trials(1,1,1,1));
control_left_incorrectc1= squeeze(sorted_trials(1,1,1,4));
control_left_incorrectc2= squeeze(sorted_trials(1,1,1,3));
control_left_incorrectc3= squeeze(sorted_trials(1,1,1,2));
control_left_incorrectc4= squeeze(sorted_trials(1,1,1,1));

led_right_correctc1 = squeeze(sorted_trials(2,2,1,4,:));
led_right_correctc2= squeeze(sorted_trials(2,2,1,3,:));
led_right_correctc3= squeeze(sorted_trials(2,2,1,2,:));
led_right_correctc4= squeeze(sorted_trials(2,2,1,1,:));
led_left_correctc1= squeeze(sorted_trials(2,2,2,4,:));
led_left_correctc2= squeeze(sorted_trials(2,2,2,3,:));                
led_left_correctc3= squeeze(sorted_trials(2,2,2,2,:));  
led_left_correctc4= squeeze(sorted_trials(2,2,2,1,:));
led_right_incorrectc1= squeeze(sorted_trials(2,1,1,4));
led_right_incorrectc2= squeeze(sorted_trials(2,1,1,3));
led_right_incorrectc3= squeeze(sorted_trials(2,1,1,2));
led_right_incorrectc4= squeeze(sorted_trials(2,1,1,1));
led_left_incorrectc1= squeeze(sorted_trials(2,1,1,4));
led_left_incorrectc2= squeeze(sorted_trials(2,1,1,3));
led_left_incorrectc3= squeeze(sorted_trials(2,1,1,2));
led_left_incorrectc4= squeeze(sorted_trials(2,1,1,1));

%%

control_idx=ones(1,length(tLeftTrial));
led_trials=find(block2);
control_idx(led_trials)=0;
led_idx=block2;
left_controls=intersect(find(control_idx), find(tLeftTrial));
right_controls=intersect(find(control_idx), find(~tLeftTrial));
left_led=intersect(find(led_idx), find(tLeftTrial));
right_led=intersect(find(led_idx), find(~tLeftTrial));
left_control_idx=zeros(1,length(tLeftTrial));
right_control_idx=zeros(1,length(tLeftTrial));
left_led_idx=zeros(1,length(tLeftTrial));
right_led_idx=zeros(1,length(tLeftTrial));
left_control_idx(left_controls)=1;
right_control_idx(right_controls)=1;
left_led_idx(left_led)=1;
right_led_idx(right_led)=1;

%%
control_c1_rightcorrecttrials_idx=nan(1,numtrials);
control_c2_rightcorrecttrials_idx=nan(1,numtrials);
control_c3_rightcorrecttrials_idx=nan(1,numtrials);
control_c4_rightcorrecttrials_idx=nan(1,numtrials);
control_c1_leftcorrecttrials_idx=nan(1,numtrials);
control_c2_leftcorrecttrials_idx=nan(1,numtrials);
control_c3_leftcorrecttrials_idx=nan(1,numtrials);
control_c4_leftcorrecttrials_idx=nan(1,numtrials);
control_c1_rightincorrecttrials_idx=nan(1,numtrials);
control_c2_rightincorrecttrials_idx=nan(1,numtrials);
control_c3_rightincorrecttrials_idx=nan(1,numtrials);
control_c4_rightincorrecttrials_idx=nan(1,numtrials);
control_c1_leftincorrecttrials_idx=nan(1,numtrials);
control_c2_leftincorrecttrials_idx=nan(1,numtrials);
control_c3_leftincorrecttrials_idx=nan(1,numtrials);
control_c4_leftincorrecttrials_idx=nan(1,numtrials);

%%
addpath('\\crash.dhe.duke.edu\rohan\GlickHull_local_analysis_code\ImagingCode-Glickfeld-Hull')
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg);
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity;
rdt=-1*ldt;
control_c1_right_idx=zeros(1,length(tLeftTrial)); %preallocate contrast+side sorted indices
control_c1_left_idx=zeros(1,length(tLeftTrial));
control_c2_right_idx=zeros(1,length(tLeftTrial));
control_c2_left_idx=zeros(1,length(tLeftTrial));
control_c3_right_idx=zeros(1,length(tLeftTrial));
control_c3_left_idx=zeros(1,length(tLeftTrial));
control_c4_right_idx=zeros(1,length(tLeftTrial));
control_c4_left_idx=zeros(1,length(tLeftTrial));

control_c1_righttrials=intersect(find(contrastratio==c1), find(right_control_idx));
control_c1_lefttrials=intersect(find(contrastratio==c1), find(left_control_idx));
control_c1_right_idx(control_c1_righttrials)=1;
control_c1_left_idx(control_c1_lefttrials)=1;
control_c1_rightcorrecttrials=intersect(find(control_c1_right_idx), find(SIx));
control_c1_leftcorrecttrials=intersect(find(control_c1_left_idx), find(SIx));
control_c1_rightincorrecttrials=intersect(find(control_c1_right_idx), find(FIx));
control_c1_leftincorrecttrials=intersect(find(control_c1_left_idx), find(FIx));

control_c2_righttrials=intersect(find(contrastratio==c2), find(right_control_idx));
control_c2_lefttrials=intersect(find(contrastratio==c2),find(left_control_idx));
control_c2_right_idx(control_c2_righttrials)=1;
control_c2_left_idx(control_c2_lefttrials)=1;
control_c2_rightcorrecttrials=intersect(find(control_c2_right_idx), find(SIx));
control_c2_leftcorrecttrials=intersect(find(control_c2_left_idx), find(SIx));
control_c2_rightincorrecttrials=intersect(find(control_c2_right_idx), find(FIx));
control_c2_leftincorrecttrials=intersect(find(control_c2_left_idx), find(FIx));

control_c3_righttrials=intersect(find(contrastratio==c3), find(right_control_idx));
control_c3_lefttrials=intersect(find(contrastratio==c3), find(left_control_idx));
control_c3_right_idx(control_c3_righttrials)=1;
control_c3_left_idx(control_c3_lefttrials)=1;
control_c3_rightcorrecttrials=intersect(find(control_c3_right_idx), find(SIx));
control_c3_leftcorrecttrials=intersect(find(control_c3_left_idx), find(SIx));
control_c3_rightincorrecttrials=intersect(find(control_c3_right_idx), find(FIx));
control_c3_leftincorrecttrials=intersect(find(control_c3_left_idx), find(FIx));

control_c4_righttrials=intersect(find(contrastratio==c4), find(right_control_idx));
control_c4_lefttrials=intersect(find(contrastratio==c4), find(left_control_idx));
control_c4_right_idx(control_c4_righttrials)=1;
control_c4_left_idx(control_c4_lefttrials)=1;
control_c4_rightcorrecttrials=intersect(find(control_c4_right_idx), find(SIx));
control_c4_leftcorrecttrials=intersect(find(control_c4_left_idx), find(SIx));
control_c4_rightincorrecttrials=intersect(find(control_c4_right_idx), find(FIx));
control_c4_leftincorrecttrials=intersect(find(control_c4_left_idx), find(FIx));
%%
for i=1:length(control_c1_rightcorrecttrials)
    control_c1_rightcorrecttrials_idx(i)=control_c1_rightcorrecttrials(i);
end
for i=1:length(control_c2_rightcorrecttrials)
    control_c2_rightcorrecttrials_idx(i)=control_c2_rightcorrecttrials(i);
end
for i=1:length(control_c3_rightcorrecttrials)
    control_c3_rightcorrecttrials_idx(i)=control_c3_rightcorrecttrials(i);
end
for i=1:length(control_c4_rightcorrecttrials)
    control_c4_rightcorrecttrials_idx(i)=control_c4_rightcorrecttrials(i);
end
for i=1:length(control_c1_leftcorrecttrials)
    control_c1_leftcorrecttrials_idx(i)=control_c1_leftcorrecttrials(i);
end
for i=1:length(control_c2_leftcorrecttrials)
    control_c2_leftcorrecttrials_idx(i)=control_c2_leftcorrecttrials(i);
end
for i=1:length(control_c3_leftcorrecttrials)
    control_c3_leftcorrecttrials_idx(i)=control_c3_leftcorrecttrials(i);
end
for i=1:length(control_c4_leftcorrecttrials)
    control_c4_leftcorrecttrials_idx(i)=control_c4_leftcorrecttrials(i);
end
for i=1:length(control_c1_rightincorrecttrials)
    control_c1_rightincorrecttrials_idx(i)=control_c1_rightincorrecttrials(i);
end
for i=1:length(control_c2_rightincorrecttrials)
    control_c2_rightincorrecttrials_idx(i)=control_c2_rightincorrecttrials(i);
end
for i=1:length(control_c3_rightincorrecttrials)
    control_c3_rightincorrecttrials_idx(i)=control_c3_rightincorrecttrials(i);
end
for i=1:length(control_c4_rightincorrecttrials)
    control_c4_rightincorrecttrials_idx(i)=control_c4_rightincorrecttrials(i);
end


%% create matrix to check for empty criteria sorted trials

empty_check_mat=[control_c1_rightcorrecttrials_idx;
    control_c2_rightcorrecttrials_idx;
    control_c3_rightcorrecttrials_idx;
    control_c4_rightcorrecttrials_idx;
    control_c1_leftcorrecttrials_idx;
    control_c2_leftcorrecttrials_idx;
    control_c3_leftcorrecttrials_idx;
    control_c4_leftcorrecttrials_idx;
    control_c1_rightincorrecttrials_idx;
    control_c2_rightincorrecttrials_idx;
    control_c3_rightincorrecttrials_idx;
    control_c4_rightincorrecttrials_idx;
    control_c1_leftincorrecttrials_idx;
    control_c2_leftincorrecttrials_idx;
    control_c3_leftincorrecttrials_idx;
    control_c4_leftincorrecttrials_idx];
        
%% control_c1_right_plot_qvals=qVals_final(:,control_c1_rightcorrecttrials)-qVals_final(8000,control_c1_rightcorrecttrials);
control_c1_right_correct_qvals=(-1)*(qVals_final(:,control_c1_rightcorrecttrials)-qVals_final(8000,control_c1_rightcorrecttrials));
control_c1_left_correct_qvals=qVals_final(:,control_c1_leftcorrecttrials)-qVals_final(8000,control_c1_leftcorrecttrials);
control_c2_right_correct_qvals=(-1)*(qVals_final(:,control_c2_rightcorrecttrials)-qVals_final(8000,control_c2_rightcorrecttrials));
control_c2_left_correct_qvals=qVals_final(:,control_c2_leftcorrecttrials)-qVals_final(8000,control_c2_leftcorrecttrials);
control_c3_right_correct_qvals=(-1)*(qVals_final(:,control_c3_rightcorrecttrials)-qVals_final(8000,control_c3_rightcorrecttrials));
control_c3_left_correct_qvals=qVals_final(:,control_c3_leftcorrecttrials)-qVals_final(8000,control_c3_leftcorrecttrials);
control_c4_right_correct_qvals=(-1)*(qVals_final(:,control_c4_rightcorrecttrials)-qVals_final(8000,control_c4_rightcorrecttrials));
control_c4_left_correct_qvals=qVals_final(:,control_c4_leftcorrecttrials)-qVals_final(8000,control_c4_leftcorrecttrials);

control_c1_left_incorrect_qvals=(-1)*(qVals_final(:,control_c1_leftincorrecttrials)-qVals_final(8000,control_c1_leftincorrecttrials));
control_c2_right_incorrect_qvals=qVals_final(:,control_c2_rightincorrecttrials)-qVals_final(8000,control_c2_rightincorrecttrials);
control_c2_left_incorrect_qvals=(-1)*(qVals_final(:,control_c2_leftincorrecttrials)-qVals_final(8000,control_c2_leftincorrecttrials));
control_c3_right_incorrect_qvals=qVals_final(:,control_c3_rightincorrecttrials)-qVals_final(8000,control_c3_rightincorrecttrials);
control_c3_left_incorrect_qvals=(-1)*(qVals_final(:,control_c3_leftincorrecttrials)-qVals_final(8000,control_c3_leftincorrecttrials));
control_c4_right_incorrect_qvals=qVals_final(:,control_c4_rightincorrecttrials)-qVals_final(8000,control_c4_rightincorrecttrials);
control_c4_left_incorrect_qvals=(-1)*(qVals_final(:,control_c4_leftincorrecttrials)-qVals_final(8000,control_c4_leftincorrecttrials));
%%
if sum(block2)~=0
    
    led_c1_right_correct_qvals=(-1)*(qVals_final(:,led_c1_rightcorrecttrials)-qVals_final(8000,led_c1_rightcorrecttrials));
    led_c1_left_correct_qvals=qVals_final(:,led_c1_leftcorrecttrials)-qVals_final(8000,led_c1_leftcorrecttrials);
    led_c2_right_correct_qvals=(-1)*(qVals_final(:,led_c2_rightcorrecttrials)-qVals_final(8000,led_c2_rightcorrecttrials));
    led_c2_left_correct_qvals=qVals_final(:,led_c2_leftcorrecttrials)-qVals_final(8000,led_c2_leftcorrecttrials);
    led_c3_right_correct_qvals=(-1)*(qVals_final(:,led_c3_rightcorrecttrials)-qVals_final(8000,led_c3_rightcorrecttrials));
    led_c3_left_correct_qvals=qVals_final(:,led_c3_leftcorrecttrials)-qVals_final(8000,led_c3_leftcorrecttrials);
    led_c4_right_correct_qvals=(-1)*(qVals_final(:,led_c4_rightcorrecttrials)-qVals_final(8000,led_c4_rightcorrecttrials));
    led_c4_left_correct_qvals=qVals_final(:,led_c4_leftcorrecttrials)-qVals_final(8000,led_c4_leftcorrecttrials);

    led_c1_left_incorrect_qvals=(-1)*(qVals_final(:,led_c1_leftincorrecttrials)-qVals_final(8000,led_c1_leftincorrecttrials));
    led_c2_right_incorrect_qvals=qVals_final(:,led_c2_rightincorrecttrials)-qVals_final(8000,led_c2_rightincorrecttrials);
    led_c2_left_incorrect_qvals=(-1)*(qVals_final(:,led_c2_leftincorrecttrials)-qVals_final(8000,led_c2_leftincorrecttrials));
    led_c3_right_incorrect_qvals=qVals_final(:,led_c3_rightincorrecttrials)-qVals_final(8000,led_c3_rightincorrecttrials);
    led_c3_left_incorrect_qvals=(-1)*(qVals_final(:,led_c3_leftincorrecttrials)-qVals_final(8000,led_c3_leftincorrecttrials));
    led_c4_right_incorrect_qvals=qVals_final(:,led_c4_rightincorrecttrials)-qVals_final(8000,lead_c4_rightincorrecttrials);
    led_c4_left_incorrect_qvals=(-1)*(qVals_final(:,led_c4_leftincorrecttrials)-qVals_final(8000,led_c4_leftincorrecttrials));
end


%%
for i=1:length(control_c1_right_correct_qvals)
control_c1_right_correct_mean_qval(i)=nanmean(control_c1_right_correct_qvals(i));
end
for i=1:length(control_c1_left_correct_qvals)
control_c1_left_correct_mean_qval(i)=nanmean(control_c1_left_correct_qvals(i));
end
for i=1:length(control_c2_right_correct_qvals)
control_c2_right_correct_mean_qval(i)=nanmean(control_c2_right_correct_qvals(i));
end
for i=1:length(control_c2_left_correct_qvals)
control_c2_left_correct_mean_qval(i)=nanmean(control_c2_left_correct_qvals(i));
end
for i=1:length(control_c3_right_correctt_qvals)
control_c3_right_correct_mean_qval(i)=nanmean(control_c3_right_correct_qvals(i));
end
for i=1:length(control_c3_left_correct_qvals)
control_c3_left_correct_mean_qval(i)=nanmean(control_c3_left_correct_qvals(i));
end
for i=1:length(control_c4_right_correct_qvals)
control_c4_right_correct_mean_qval(i)=nanmean(control_c4_right_correct_qvals(i));
end
for i=1:length(control_c4_left_correct_qvals)
control_c4_left_correct_mean_qval(i)=nanmean(control_c4_left_correct_qvals(i));
end
% repeat for incorrect trials below 
for i=1:length(control_c1_right_incorrect_qvals)
control_c1_right_incorrect_mean_qval(i)=nanmean(control_c1_right_incorrect_qvals(i));
end
for i=1:length(control_c1_left_incorrect_qvals)
control_c1_left_incorrect_mean_qval(i)=nanmean(control_c1_left_incorrect_qvals(i));
end
for i=1:length(control_c2_right_incorrect_qvals)
control_c2_right_incorrect_mean_qval(i)=nanmean(control_c2_right_incorrect_qvals(i));
end
for i=1:length(control_c2_left_incorrect_qvals)
control_c2_left_incorrect_mean_qval(i)=nanmean(control_c2_left_incorrect_qvals(i));
end
for i=1:length(control_c3_right_incorrectt_qvals)
control_c3_right_incorrect_mean_qval(i)=nanmean(control_c3_right_incorrect_qvals(i));
end
for i=1:length(control_c3_left_incorrect_qvals)
control_c3_left_incorrect_mean_qval(i)=nanmean(control_c3_left_incorrect_qvals(i));
end
for i=1:length(control_c4_right_incorrect_qvals)
control_c4_right_incorrect_mean_qval(i)=nanmean(control_c4_right_incorrect_qvals(i));
end
for i=1:length(control_c4_left_incorrect_qvals)
control_c4_left_incorrect_mean_qval(i)=nanmean(control_c4_left_incorrect_qvals(i));
end
%%
if sum(block2)~=0
    for i=1:length(led_c1_right_correct_qvals)
    led_c1_right_correct_mean_qval(i)=nanmean(led_c1_right_correct_qvals(i));
    end
    for i=1:length(led_c1_left_correct_qvals)
    led_c1_left_correct_mean_qval(i)=nanmean(led_c1_left_correct_qvals(i));
    end
    for i=1:length(led_c2_right_correct_qvals)
    led_c2_right_correct_mean_qval(i)=nanmean(led_c2_right_correct_qvals(i));
    end
    for i=1:length(led_c2_left_correct_qvals)
    led_c2_left_correct_mean_qval(i)=nanmean(led_c2_left_correct_qvals(i));
    end
    for i=1:length(led_c3_right_correct_qvals)
    led_c3_right_correct_mean_qval(i)=nanmean(led_c3_right_correct_qvals(i));
    end
    for i=1:length(led_c3_left_correct_qvals)
    led_c3_left_correct_mean_qval(i)=nanmean(led_c3_left_correct_qvals(i));
    end
    for i=1:length(led_c4_right_correct_qvals)
    led_c4_right_correct_mean_qval(i)=nanmean(led_c4_right_correct_qvals(i));
    end
    for i=1:length(led_c4_left_correct_qvals)
    led_c4_left_correct_mean_qval(i)=nanmean(led_c4_left_correct_qvals(i));
    end
    % repeat for incorrect trials below 
    for i=1:length(led_c1_right_incorrect_qvals)
    led_c1_right_incorrect_mean_qval(i)=nanmean(led_c1_right_incorrect_qvals(i));
    end
    for i=1:length(led_c1_left_incorrect_qvals)
    led_c1_left_incorrect_mean_qval(i)=nanmean(led_c1_left_incorrect_qvals(i));
    end
    for i=1:length(led_c2_right_incorrect_qvals)
    led_c2_right_incorrect_mean_qval(i)=nanmean(led_c2_right_incorrect_qvals(i));
    end
    for i=1:length(led_c2_left_incorrect_qvals)
    led_c2_left_incorrect_mean_qval(i)=nanmean(led_c2_left_incorrect_qvals(i));
    end
    for i=1:length(led_c3_right_incorrect_qvals)
    led_c3_right_incorrect_mean_qval(i)=nanmean(led_c3_right_incorrect_qvals(i));
    end
    for i=1:length(led_c3_left_incorrect_qvals)
    led_c3_left_incorrect_mean_qval(i)=nanmean(led_c3_left_incorrect_qvals(i));
    end
    for i=1:length(led_c4_right_incorrect_qvals)
    led_c4_right_incorrect_mean_qval(i)=nanmean(led_c4_right_incorrect_qvals(i));
    end
    for i=1:length(led_c4_left_incorrect_qvals)
    led_c4_left_incorrect_mean_qval(i)=nanmean(led_c4_left_incorrect_qvals(i));
    end
end
control_mean_qvals_mat_top=[control_c1_left_correct_mean_qval;
    control_c1_right_incorrect_mean_qval;
    control_c2_left_correct_mean_qval;
    control_c2_right_incorrect_mean_qval;
    control_c3_left_correct_mean_qval;
    control_c3_right_incorrect_mean_qval;
    control_c4_left_correct_mean_qval;
    control_c4_right_incorrect_mean_qval;];
control_mean_qvals_mat_bottom=[control_c1_left_incorrect_mean_qval;
    control_c1_right_correct_mean_qval;
    control_c2_left_incorrect_mean_qval;
    control_c2_right_correct_mean_qval;
    control_c3_left_incorrect_mean_qval;
    control_c3_right_correct_mean_qval;
    control_c4_left_incorrect_mean_qval;
    control_c4_right_correct_mean_qval;];
if sum(block2)~=0
    led_mean_qvals_mat_top=[led_c1_left_correct_mean_qval;
        led_c1_right_incorrect_mean_qval;
        led_c2_left_correct_mean_qval;
        led_c2_right_incorrect_mean_qval;
        led_c3_left_correct_mean_qval;
        led_c3_right_incorrect_mean_qval;
        led_c4_left_correct_mean_qval;
        led_c4_right_incorrect_mean_point_qval;];
    
    led_mean_qvals_mat_bottom=[led_c1_left_incorrect_mean_qval;
        led_c1_right_correct_mean_qval;
        led_c2_left_incorrect_mean_qval;
        led_c2_right_correct_mean_qval;
        led_c3_left_incorrect_mean_qval;
        led_c3_right_correct_mean_qval;
        led_c4_left_incorrect_mean_qval;
        led_c4_right_correct_mean_point_qval;];
end

plot(qTimes_final, control_mean_qvals_mat_top)
