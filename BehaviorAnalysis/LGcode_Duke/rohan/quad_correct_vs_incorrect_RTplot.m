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
control_idx=ones(1,length(tLeftTrial));
led_trials=find(block2);
control_idx(led_trials)=0;
led_idx=block2;
left_controls=intersect(find(control_idx), find(tLeftTrial));
right_controls=intersect(find(control_idx), find(~tLeftTrial));
% left_led=intersect(find(led_idx), find(tLeftTrial));
% right_led=intersect(find(led_idx), find(~tLeftTrial));
% left_control_idx=zeros(1,length(tLeftTrial));
% right_control_idx=zeros(1,length(tLeftTrial));
% left_led_idx=zeros(1,length(tLeftTrial));
% right_led_idx=zeros(1,length(tLeftTrial));
% left_control_idx(left_controls)=1;
% right_control_idx(right_controls)=1;
% left_led_idx(left_led)=1;
% right_led_idx(right_led)=1;
%%
addpath('\\crash.dhe.duke.edu\rohan\GlickHull_local_analysis_code\ImagingCode-Glickfeld-Hull')
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg);
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity;
rdt=-1*ldt;
c1_correct_idx=zeros(1,length(tLeftTrial)); %preallocate contrast+side sorted indices
c1_correct_idx=zeros(1,length(tLeftTrial));
c2_correct_idx=zeros(1,length(tLeftTrial));
c2_correct_idx=zeros(1,length(tLeftTrial));
c3_correct_idx=zeros(1,length(tLeftTrial));
c3_correct_idx=zeros(1,length(tLeftTrial));
c4_correct_idx=zeros(1,length(tLeftTrial));
c4_correct_idx=zeros(1,length(tLeftTrial));

c1_incorrect_idx=zeros(1,length(tLeftTrial)); %preallocate contrast+side sorted indices
c1_incorrect_idx=zeros(1,length(tLeftTrial));
c2_incorrect_idx=zeros(1,length(tLeftTrial));
c2_incorrect_idx=zeros(1,length(tLeftTrial));
c3_incorrect_idx=zeros(1,length(tLeftTrial));
c3_incorrect_idx=zeros(1,length(tLeftTrial));
c4_incorrect_idx=zeros(1,length(tLeftTrial));
c4_incorrect_idx=zeros(1,length(tLeftTrial));

%%
c1_correct=intersect(find(contrastratio==c1), find(SIx));
c2_correct=intersect(find(contrastratio==c2), find(SIx));
c3_correct=intersect(find(contrastratio==c3), find(SIx));
c4_correct=intersect(find(contrastratio==c4), find(SIx));

c1_incorrect=intersect(find(contrastratio==c1), find(FIx));
c2_incorrect=intersect(find(contrastratio==c2), find(FIx));
c3_incorrect=intersect(find(contrastratio==c3), find(FIx));
c4_incorrect=intersect(find(contrastratio==c4), find(FIx));

    
c1_correct_idx(c1_correct)=1;
c2_correct_idx(c2_correct)=1;
c3_correct_idx(c3_correct)=1;
c4_correct_idx(c4_correct)=1;

c1_incorrect_idx(c1_incorrect)=1;
c2_incorrect_idx(c2_incorrect)=1;
c3_incorrect_idx(c3_incorrect)=1;
c4_incorrect_idx(c4_incorrect)=1;



%%
control_c1_correct_plot_qvals=qVals_final(:,c1_correct)-qVals_final(8000,c1_correct);
control_c2_correct_plot_qvals=qVals_final(:,c2_correct)-qVals_final(8000,c2_correct);
control_c3_correct_plot_qvals=qVals_final(:,c3_correct)-qVals_final(8000,c3_correct);
control_c4_correct_plot_qvals=qVals_final(:,c4_correct)-qVals_final(8000,c4_correct);

control_c1_incorrect_plot_qvals=qVals_final(:,c1_incorrect)-qVals_final(8000,c1_incorrect);
control_c2_incorrect_plot_qvals=qVals_final(:,c2_incorrect)-qVals_final(8000,c2_incorrect);
control_c3_incorrect_plot_qvals=qVals_final(:,c3_incorrect)-qVals_final(8000,c3_incorrect);
control_c4_incorrect_plot_qvals=qVals_final(:,c4_incorrect)-qVals_final(8000,c4_incorrect);

for i=1:length(control_c1_correct_plot_qvals)
control_c1_correct_mean_qval(i)=nanmean(control_c1_correct_plot_qvals(i));
end
for i=1:length(control_c2_correct_plot_qvals)
control_c2_correct_mean_qval(i)=nanmean(control_c2_correct_plot_qvals(i));
end
for i=1:length(control_c3_correct_plot_qvals)
control_c3_correct_mean_qval(i)=nanmean(control_c3_correct_plot_qvals(i));
end
for i=1:length(control_c4_correct_plot_qvals)
control_c4_correct_mean_qval(i)=nanmean(control_c4_correct_plot_qvals(i));
end

for i=1:length(control_c1_incorrect_plot_qvals)
control_c1_incorrect_mean_qval(i)=nanmean(control_c1_incorrect_plot_qvals(i));
end
for i=1:length(control_c2_incorrect_plot_qvals)
control_c2_incorrect_mean_qval(i)=nanmean(control_c2_incorrect_plot_qvals(i));
end
for i=1:length(control_c3_incorrect_plot_qvals)
control_c3_incorrect_mean_qval(i)=nanmean(control_c3_incorrect_plot_qvals(i));
end
for i=1:length(control_c4_incorrect_plot_qvals)
control_c4_incorrect_mean_qval(i)=nanmean(control_c4_incorrect_plot_qvals(i));
end

mean_qvals_mat=[control_c1_correct_mean_qval;
    control_c2_correct_mean_qval;
    control_c3_correct_mean_qval;
    control_c4_correct_mean_qval;
    control_c1_incorrect_mean_qval;
    control_c2_incorrect_mean_qval;
    control_c3_incorrect_mean_qval;
    control_c4_incorrect_mean_qval;];

colors=brewermap(10, 'RdBu');
colors2=brewermap(10, '*RdBu');

control_c1_correct_postTFT=control_c1_correct_mean_qval-control_c1_correct_mean_qval(8800);
absval_control_c1_correct_postTFT=abs(control_c1_correct_postTFT);
control_c1_correct_postTFT_thresh=find(absval_control_c1_correct_postTFT(8800:18000)>=ldt)+8800;
control_c1_correct_postTFT_dec=min(control_c1_correct_postTFT_thresh)-8000;

control_c2_correct_postTFT=control_c2_correct_mean_qval-control_c2_correct_mean_qval(8800);
absval_control_c2_correct_postTFT=abs(control_c2_correct_postTFT);
control_c2_correct_postTFT_thresh=find(absval_control_c2_correct_postTFT(8800:18000)>=ldt)+8800;
control_c2_correct_postTFT_dec=min(control_c2_correct_postTFT_thresh)-8000;

control_c3_correct_postTFT=control_c3_correct_mean_qval-control_c3_correct_mean_qval(8800);
absval_control_c3_correct_postTFT=abs(control_c3_correct_postTFT);
control_c3_correct_postTFT_thresh=find(absval_control_c3_correct_postTFT(8800:18000)>=ldt)+8800;
control_c3_correct_postTFT_dec=min(control_c3_correct_postTFT_thresh)-8000;

control_c4_correct_postTFT=control_c4_correct_mean_qval-control_c4_correct_mean_qval(8800);
absval_control_c4_correct_postTFT=abs(control_c4_correct_postTFT);
control_c4_correct_postTFT_thresh=find(absval_control_c4_correct_postTFT(8800:18000)>=ldt)+8800;
control_c4_correct_postTFT_dec=min(control_c4_correct_postTFT_thresh)-8000;

control_c1_incorrect_postTFT=control_c1_incorrect_mean_qval-control_c1_incorrect_mean_qval(8800);
absval_control_c1_incorrect_postTFT=abs(control_c1_incorrect_postTFT);
control_c1_incorrect_postTFT_thresh=find(absval_control_c1_incorrect_postTFT(8800:18000)>=ldt)+8800;
control_c1_incorrect_postTFT_dec=min(control_c1_incorrect_postTFT_thresh)-8000;

control_c2_incorrect_postTFT=control_c2_incorrect_mean_qval-control_c2_incorrect_mean_qval(8800);
absval_control_c2_incorrect_postTFT=abs(control_c2_incorrect_postTFT);
control_c2_incorrect_postTFT_thresh=find(absval_control_c2_incorrect_postTFT(8800:18000)>=ldt)+8800;
control_c2_incorrect_postTFT_dec=min(control_c2_incorrect_postTFT_thresh)-8000;

control_c3_incorrect_postTFT=control_c3_incorrect_mean_qval-control_c3_incorrect_mean_qval(8800);
absval_control_c3_incorrect_postTFT=abs(control_c3_incorrect_postTFT);
control_c3_incorrect_postTFT_thresh=find(absval_control_c3_incorrect_postTFT(8800:18000)>=ldt)+8800;
control_c3_incorrect_postTFT_dec=min(control_c3_incorrect_postTFT_thresh)-8000;

control_c4_incorrect_postTFT=control_c4_incorrect_mean_qval-control_c4_incorrect_mean_qval(8800);
absval_control_c4_incorrect_postTFT=abs(control_c4_incorrect_postTFT);
control_c4_incorrect_postTFT_thresh=find(absval_control_c4_incorrect_postTFT(8800:18000)>=ldt)+8800;
control_c4_incorrect_postTFT_dec=min(control_c4_incorrect_postTFT_thresh)-8000;

%%
for i=1:4
        tempvar=find(mean_qvals_mat(i,8000:18000)<= rdt| mean_qvals_mat(i,8000:18000)>=ldt);
         control_correct_dectime(i)=min(tempvar);
end
for i=5:8
    tempvar=find(mean_qvals_mat(i,8000:18000)>= ldt| mean_qvals_mat(i,8000:18000)<=rdt);
    control_incorrect_dectime(i-4)=min(tempvar);
end

if input.doSizeDiscrim==1
    taskstr='size';
end
if input.doContrastDiscrim==1
    taskstr='contrast';
end
num_trials=length(tLeftTrial);
numignores=(length(tLeftTrial)-(sum(FIx)+sum(SIx)));
num_includedtrials=num_trials-numignores;
numcorrect=sum(SIx);
numincorrect=sum(FIx);
%%
figure;
hold on
ylim([0 1000])
pc=plot(uncontrasts, control_correct_dectime,'-greeno','LineWidth',1.5, 'MarkerEdgeColor','k', 'MarkerFaceColor',colors2(4,:), 'MarkerSize',5)
pi=plot(uncontrasts, control_incorrect_dectime, '-redo', 'LineWidth', 1.5, 'MarkerEdgeColor','k','MarkerFaceColor', 'red', 'MarkerSize', 5)
pc.Color=colors2(4,:);
xstr=sprintf('Target/Distractor %s ratio', taskstr);
xlabel(xstr); ylabel('Reaction time (ms)')
titlestr=sprintf('i435 Average reaction time for correct and incorrect trials by %s ratio in example session %s \nN=%d trials. %d correct, %d incorrect', taskstr, savedate(1:6), num_includedtrials, numcorrect, numincorrect);
title(titlestr);
legend('Correct trials', 'Incorrect trials')
%%
colors=brewermap(13, '*Reds');
colors2=brewermap(13,'*Greens');
figure;
hold on
patch([-2000 0 0 -2000],[-1000 -1000 1000 1000],[0.2 0.2 0.2])
patch([ 0 800 800 0],[-1000 -1000 1000 1000],[0.7 0.7 0.8])
set(gca,'children',flipud(get(gca,'children')))
rdtline=line([-8000 10000], [rdt rdt], 'Color', 'yellow','LineWidth', 3,'LineStyle','--')
ldtline=line([-8000 10000], [ldt ldt], 'Color', 'green', 'LineWidth', 3, 'LineStyle','--')
xlim([-500 1500])
ylim([-200 200])
p1=plot(qTimes_final, abs(mean_qvals_mat(1,:)),'LineWidth', 1.7)
p2=plot(qTimes_final, abs(mean_qvals_mat(2,:)),'LineWidth', 1.7)
p3=plot(qTimes_final, abs(mean_qvals_mat(3,:)),'LineWidth', 1.7)
p4=plot(qTimes_final, abs(mean_qvals_mat(4,:)),'LineWidth', 1.7)
p5=plot(qTimes_final, -1*abs(mean_qvals_mat(5,:)),'LineWidth', 1.7)
p6=plot(qTimes_final, -1*abs(mean_qvals_mat(6,:)),'LineWidth', 1.7)
p7=plot(qTimes_final, -1*abs(mean_qvals_mat(7,:)),'LineWidth', 1.7)
p8=plot(qTimes_final, -1*abs(mean_qvals_mat(8,:)),'LineWidth', 1.7)
xlabel('Trial time aligned to stimulus on (ms)'); ylabel('Absolute value of Quadrature Position')
titlestr=sprintf('i435 Average Quadrature Motion Paths by correct vs. incorrect trials and %s ratio in example session %s', taskstr, savedate(1:6));
title(titlestr)
legend('Stationary Period', 'Too Fast time', 'incorrect decision threshold', 'Correct decision threshold',...
    'Ratio 1 correct', 'Ratio 2 correct', 'Ratio 3 correct', 'Ratio 4 correct',...
     'Ratio 1 incorrect', 'Ratio 2 incorrect', 'Ratio 3 incorrect', 'Ratio 4 incorrect')
p1.Color=colors2(2,:);
p2.Color=colors2(4,:);
p3.Color=colors2(6,:);
p4.Color=colors2(8,:);
p5.Color=colors(2,:);
p6.Color=colors(4,:);
p7.Color=colors(6,:);
p8.Color=colors(8,:);






