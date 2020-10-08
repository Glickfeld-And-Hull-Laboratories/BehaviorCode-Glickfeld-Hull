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
temp = squeeze(sorted_trials(1,2,1,4,:))';
temp = temp(find(~isnan(temp)));
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1
control_right_correctc1 = squeeze(sorted_trials(1,2,1,4,:))';
control_right_correctc2= squeeze(sorted_trials(1,2,1,3,:))';
control_right_correctc3= squeeze(sorted_trials(1,2,1,2,:))';
control_right_correctc4= squeeze(sorted_trials(1,2,1,1,:))';
control_left_correctc1= squeeze(sorted_trials(1,2,2,4,:))';
control_left_correctc2= squeeze(sorted_trials(1,2,2,3,:))';                
control_left_correctc3= squeeze(sorted_trials(1,2,2,2,:))';  
control_left_correctc4= squeeze(sorted_trials(1,2,2,1,:))';
control_right_incorrectc1= squeeze(sorted_trials(1,1,1,4))';
control_right_incorrectc2= squeeze(sorted_trials(1,1,1,3))';
control_right_incorrectc3= squeeze(sorted_trials(1,1,1,2))';
control_right_incorrectc4= squeeze(sorted_trials(1,1,1,1))';
control_left_incorrectc1= squeeze(sorted_trials(1,1,1,4))';
control_left_incorrectc2= squeeze(sorted_trials(1,1,1,3))';
control_left_incorrectc3= squeeze(sorted_trials(1,1,1,2))';
control_left_incorrectc4= squeeze(sorted_trials(1,1,1,1))';

led_right_correctc1 = squeeze(sorted_trials(2,2,1,4,:))';
led_right_correctc2= squeeze(sorted_trials(2,2,1,3,:))';
led_right_correctc3= squeeze(sorted_trials(2,2,1,2,:))';
led_right_correctc4= squeeze(sorted_trials(2,2,1,1,:))';
led_left_correctc1= squeeze(sorted_trials(2,2,2,4,:))';
led_left_correctc2= squeeze(sorted_trials(2,2,2,3,:))';                
led_left_correctc3= squeeze(sorted_trials(2,2,2,2,:))';  
led_left_correctc4= squeeze(sorted_trials(2,2,2,1,:))';
led_right_incorrectc1= squeeze(sorted_trials(2,1,1,4))';
led_right_incorrectc2= squeeze(sorted_trials(2,1,1,3))';
led_right_incorrectc3= squeeze(sorted_trials(2,1,1,2))';
led_right_incorrectc4= squeeze(sorted_trials(2,1,1,1))';
led_left_incorrectc1= squeeze(sorted_trials(2,1,1,4))';
led_left_incorrectc2= squeeze(sorted_trials(2,1,1,3))';
led_left_incorrectc3= squeeze(sorted_trials(2,1,1,2))';
led_left_incorrectc4= squeeze(sorted_trials(2,1,1,1))';
%%
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg);
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity;
rdt=-1*ldt;
%%
control_c1_correct_right_qvals=qVals_final(:,find(~isnan(control_right_correctc1)))-qVals_final(8000,find(~isnan(control_right_correctc1)));
control_c2_correct_right_qvals=qVals_final(:,find(~isnan(control_right_correctc2)))-qVals_final(8000,find(~isnan(control_right_correctc2)));
control_c3_correct_right_qvals=qVals_final(:,find(~isnan(control_right_correctc3)))-qVals_final(8000,find(~isnan(control_right_correctc3)));
control_c4_correct_right_qvals=qVals_final(:,find(~isnan(control_right_correctc4)))-qVals_final(8000,find(~isnan(control_right_correctc4)));

control_c1_incorrect_right_qvals=qVals_final(:,find(~isnan(control_right_incorrectc1)))-qVals_final(8000,find(~isnan(control_right_incorrectc1)));
control_c2_incorrect_right_qvals=qVals_final(:,find(~isnan(control_right_incorrectc2)))-qVals_final(8000,find(~isnan(control_right_incorrectc2)));
control_c3_incorrect_right_qvals=qVals_final(:,find(~isnan(control_right_incorrectc3)))-qVals_final(8000,find(~isnan(control_right_incorrectc3)));
control_c4_incorrect_right_qvals=qVals_final(:,find(~isnan(control_right_incorrectc4)))-qVals_final(8000,find(~isnan(control_right_incorrectc4)));

control_c1_correct_left_qvals=qVals_final(:,find(~isnan(control_left_correctc1)))-qVals_final(8000,find(~isnan(control_left_correctc1)));
control_c2_correct_left_qvals=qVals_final(:,find(~isnan(control_left_correctc2)))-qVals_final(8000,find(~isnan(control_left_correctc2)));
control_c3_correct_left_qvals=qVals_final(:,find(~isnan(control_left_correctc3)))-qVals_final(8000,find(~isnan(control_left_correctc3)));
control_c4_correct_left_qvals=qVals_final(:,find(~isnan(control_left_correctc4)))-qVals_final(8000,find(~isnan(control_left_correctc4)));

control_c1_incorrect_left_qvals=qVals_final(:,find(~isnan(control_left_incorrectc1)))-qVals_final(8000,find(~isnan(control_left_incorrectc1)));
control_c2_incorrect_left_qvals=qVals_final(:,find(~isnan(control_left_incorrectc2)))-qVals_final(8000,find(~isnan(control_left_incorrectc2)));
control_c3_incorrect_left_qvals=qVals_final(:,find(~isnan(control_left_incorrectc3)))-qVals_final(8000,find(~isnan(control_left_incorrectc3)));
control_c4_incorrect_left_qvals=qVals_final(:,find(~isnan(control_left_incorrectc4)))-qVals_final(8000,find(~isnan(control_left_incorrectc4)));

led_c1_correct_right_qvals=qVals_final(:,find(~isnan(led_right_correctc1)))-qVals_final(8000,find(~isnan(led_right_correctc1)));
led_c2_correct_right_qvals=qVals_final(:,find(~isnan(led_right_correctc2)))-qVals_final(8000,find(~isnan(led_right_correctc2)));
led_c3_correct_right_qvals=qVals_final(:,find(~isnan(led_right_correctc3)))-qVals_final(8000,find(~isnan(led_right_correctc3)));
led_c4_correct_right_qvals=qVals_final(:,find(~isnan(led_right_correctc4)))-qVals_final(8000,find(~isnan(led_right_correctc4)));

led_c1_incorrect_right_qvals=qVals_final(:,find(~isnan(led_right_incorrectc1)))-qVals_final(8000,find(~isnan(led_right_incorrectc1)));
led_c2_incorrect_right_qvals=qVals_final(:,find(~isnan(led_right_incorrectc2)))-qVals_final(8000,find(~isnan(led_right_incorrectc2)));
led_c3_incorrect_right_qvals=qVals_final(:,find(~isnan(led_right_incorrectc3)))-qVals_final(8000,find(~isnan(led_right_incorrectc3)));
led_c4_incorrect_right_qvals=qVals_final(:,find(~isnan(led_right_incorrectc4)==1))-qVals_final(8000,find(~isnan(led_right_incorrectc4)));

led_c1_correct_left_qvals=qVals_final(:,find(~isnan(led_left_correctc1)))-qVals_final(8000,find(~isnan(led_left_correctc1)));
led_c2_correct_left_qvals=qVals_final(:,find(~isnan(led_left_correctc2)))-qVals_final(8000,find(~isnan(led_left_correctc2)));
led_c3_correct_left_qvals=qVals_final(:,find(~isnan(led_left_correctc3)))-qVals_final(8000,find(~isnan(led_left_correctc3)));
led_c4_correct_left_qvals=qVals_final(:,find(~isnan(led_left_correctc4)))-qVals_final(8000,find(~isnan(led_left_correctc4)));

led_c1_incorrect_left_qvals=qVals_final(:,find(~isnan(led_left_incorrectc1)))-qVals_final(8000,find(~isnan(led_left_incorrectc1)));
led_c2_incorrect_left_qvals=qVals_final(:,find(~isnan(led_left_incorrectc2)))-qVals_final(8000,find(~isnan(led_left_incorrectc2)));
led_c3_incorrect_left_qvals=qVals_final(:,find(~isnan(led_left_incorrectc3)))-qVals_final(8000,find(~isnan(led_left_incorrectc3)));
led_c4_incorrect_left_qvals=qVals_final(:,find(~isnan(led_left_incorrectc4)))-qVals_final(8000,find(~isnan(led_left_incorrectc4)));

colors=brewermap(13, '*Reds');
colors2=brewermap(13,'*Greens');

figure;
 
patch([-2000 0 0 -2000],[-1000 -1000 1000 1000],[0.2 0.2 0.2])
patch([ 0 800 800 0],[-1000 -1000 1000 1000],[0.7 0.7 0.8])
set(gca,'children',flipud(get(gca,'children')))
rdtline=line([-8000 10000], [rdt rdt], 'Color', 'yellow','LineWidth', 3,'LineStyle','--')
ldtline=line([-8000 10000], [ldt ldt], 'Color', 'green', 'LineWidth', 3, 'LineStyle','--')
ylim([-300 400])
xlim([-500 2000]) 
hold on
p1=plot(qTimes_final, control_c1_correct_right_qvals, 'Color', colors(2,:))
p2=plot(qTimes_final, control_c2_correct_right_qvals)
p3=plot(qTimes_final, control_c3_correct_right_qvals)
p4=plot(qTimes_final, control_c4_correct_right_qvals)
p5=plot(qTimes_final, control_c1_incorrect_left_qvals)
p6=plot(qTimes_final, control_c2_incorrect_left_qvals)
p7=plot(qTimes_final, control_c3_incorrect_left_qvals)
p8=plot(qTimes_final, control_c4_incorrect_left_qvals)
%%
p1.Color=colors(1,:)
p2.Color=colors2(1,:)
p3.Color=colors(2,:)
p4.Color=colors2(2,:)
p5.Color=colors(3,:)
p6.Color=colors2(2,:)
p7.Color=colors(4,:)
p8.Color=colors2(4,:)
