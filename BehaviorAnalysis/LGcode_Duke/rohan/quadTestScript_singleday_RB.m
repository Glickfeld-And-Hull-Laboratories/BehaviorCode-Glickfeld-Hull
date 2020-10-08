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
        qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]);
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

%% need to sort into control and LED
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

control_c2_righttrials=intersect(find(contrastratio==c2), find(right_control_idx));
control_c2_lefttrials=intersect(find(contrastratio==c2),find(left_control_idx));
control_c2_right_idx(control_c2_righttrials)=1;
control_c2_left_idx(control_c2_lefttrials)=1;
control_c2_rightcorrecttrials=intersect(find(control_c2_right_idx), find(SIx));
control_c2_leftcorrecttrials=intersect(find(control_c2_left_idx), find(SIx));

control_c3_righttrials=intersect(find(contrastratio==c3), find(right_control_idx));
control_c3_lefttrials=intersect(find(contrastratio==c3), find(left_control_idx));
control_c3_right_idx(control_c3_righttrials)=1;
control_c3_left_idx(control_c3_lefttrials)=1;
control_c3_rightcorrecttrials=intersect(find(control_c3_right_idx), find(SIx));
control_c3_leftcorrecttrials=intersect(find(control_c3_left_idx), find(SIx));

control_c4_righttrials=intersect(find(contrastratio==c4), find(right_control_idx));
control_c4_lefttrials=intersect(find(contrastratio==c4), find(left_control_idx));
control_c4_right_idx(control_c4_righttrials)=1;
control_c4_left_idx(control_c4_lefttrials)=1;
control_c4_rightcorrecttrials=intersect(find(control_c4_right_idx), find(SIx));
control_c4_leftcorrecttrials=intersect(find(control_c4_left_idx), find(SIx));

control_c1_right_plot_qvals=qVals_final(:,control_c1_rightcorrecttrials)-qVals_final(8000,control_c1_rightcorrecttrials);
control_c1_left_plot_qvals=qVals_final(:,control_c1_leftcorrecttrials)-qVals_final(8000,control_c1_leftcorrecttrials);
control_c2_right_plot_qvals=qVals_final(:,control_c2_rightcorrecttrials)-qVals_final(8000,control_c2_rightcorrecttrials);
control_c2_left_plot_qvals=qVals_final(:,control_c2_leftcorrecttrials)-qVals_final(8000,control_c2_leftcorrecttrials);
control_c3_right_plot_qvals=qVals_final(:,control_c3_rightcorrecttrials)-qVals_final(8000,control_c3_rightcorrecttrials);
control_c3_left_plot_qvals=qVals_final(:,control_c3_leftcorrecttrials)-qVals_final(8000,control_c3_leftcorrecttrials);
control_c4_right_plot_qvals=qVals_final(:,control_c4_rightcorrecttrials)-qVals_final(8000,control_c4_rightcorrecttrials);
control_c4_left_plot_qvals=qVals_final(:,control_c4_leftcorrecttrials)-qVals_final(8000,control_c4_leftcorrecttrials);

for i=1:length(control_c1_right_plot_qvals)
control_c1_right_mean_point_qval(i)=nanmean(control_c1_right_plot_qvals(i));
end
for i=1:length(control_c1_left_plot_qvals)
control_c1_left_mean_point_qval(i)=nanmean(control_c1_left_plot_qvals(i));
end
for i=1:length(control_c2_right_plot_qvals)
control_c2_right_mean_point_qval(i)=nanmean(control_c2_right_plot_qvals(i));
end
for i=1:length(control_c2_left_plot_qvals)
control_c2_left_mean_point_qval(i)=nanmean(control_c2_left_plot_qvals(i));
end
for i=1:length(control_c3_right_plot_qvals)
control_c3_right_mean_point_qval(i)=nanmean(control_c3_right_plot_qvals(i));
end
for i=1:length(control_c3_left_plot_qvals)
control_c3_left_mean_point_qval(i)=nanmean(control_c3_left_plot_qvals(i));
end
for i=1:length(control_c4_right_plot_qvals)
control_c4_right_mean_point_qval(i)=nanmean(control_c4_right_plot_qvals(i));
end
for i=1:length(control_c4_left_plot_qvals)
control_c4_left_mean_point_qval(i)=nanmean(control_c4_left_plot_qvals(i));
end
mean_qvals_mat=[control_c1_right_mean_point_qval;
    control_c2_right_mean_point_qval;
    control_c3_right_mean_point_qval;
    control_c4_right_mean_point_qval;
    control_c1_left_mean_point_qval;
    control_c2_left_mean_point_qval;
    control_c3_left_mean_point_qval;
    control_c4_left_mean_point_qval;];
colors=brewermap(10, 'RdBu');
colors2=brewermap(10, '*RdBu');

control_c1_rightpostTFT=control_c1_right_mean_point_qval-control_c1_right_mean_point_qval(8800);
control_c1_rightpostTFT_thresh=find(control_c1_rightpostTFT(8800:18000)<=rdt)+8800;
control_c1_rightpostTFT_dec=min(control_c1_rightpostTFT_thresh)-8000;

control_c2_rightpostTFT=control_c2_right_mean_point_qval-control_c2_right_mean_point_qval(8800);
control_c2_rightpostTFT_thresh=find(control_c2_rightpostTFT(8800:18000)<=rdt)+8800;
control_c2_rightpostTFT_dec=min(control_c2_rightpostTFT_thresh)-8000;
if isempty(control_c2_rightpostTFT_thresh)
    control_c2_rightpostTFT_dec=NaN;
end

control_c3_rightpostTFT=control_c3_right_mean_point_qval-control_c3_right_mean_point_qval(8800);
control_c3_rightpostTFT_thresh=find(control_c3_rightpostTFT(8800:18000)<=rdt)+8800;
control_c3_rightpostTFT_dec=min(control_c3_rightpostTFT_thresh)-8000;

control_c4_rightpostTFT=control_c4_right_mean_point_qval-control_c4_right_mean_point_qval(8800);
control_c4_rightpostTFT_thresh=find(control_c4_rightpostTFT(8800:18000)<=rdt)+8800;
control_c4_rightpostTFT_dec=min(control_c4_rightpostTFT_thresh)-8000;

control_allright_postTFT_dec=[control_c1_rightpostTFT_dec, control_c2_rightpostTFT_dec, control_c3_rightpostTFT_dec, control_c4_rightpostTFT_dec];
%same for left side below

control_c1_leftpostTFT=control_c1_left_mean_point_qval-control_c1_left_mean_point_qval(8800);
control_c1_leftpostTFT_thresh=find(control_c1_leftpostTFT(8800:18000)>=ldt)+8800;
control_c1_leftpostTFT_dec=min(control_c1_leftpostTFT_thresh)-8000; 

control_c2_leftpostTFT=control_c2_left_mean_point_qval-control_c2_left_mean_point_qval(8800);
control_c2_leftpostTFT_thresh=find(control_c2_leftpostTFT(8800:18000)>=ldt)+8800;
control_c2_leftpostTFT_dec=min(control_c2_leftpostTFT_thresh)-8000;

control_c3_leftpostTFT=control_c3_left_mean_point_qval-control_c3_left_mean_point_qval(8800);
control_c3_leftpostTFT_thresh=find(control_c3_leftpostTFT(8800:18000)>=ldt)+8800;
control_c3_leftpostTFT_dec=min(control_c3_leftpostTFT_thresh)-8000;

control_c4_leftpostTFT=control_c4_left_mean_point_qval-control_c4_left_mean_point_qval(8800);
control_c4_leftpostTFT_thresh=find(control_c4_leftpostTFT(8800:18000)>=ldt)+8800;
control_c4_leftpostTFT_dec=min(control_c4_leftpostTFT_thresh)-8000;

control_allleft_postTFT_dec=[control_c1_leftpostTFT_dec, control_c2_leftpostTFT_dec, control_c3_leftpostTFT_dec, control_c4_leftpostTFT_dec];

for i=1:4
    control_right_actual_dec(i)=mean_qvals_mat(i, control_allright_postTFT_dec(i)+8000); %actual quadrature positions at measured RT's
end
for i=1:4
    control_left_actual_dec(i)=mean_qvals_mat(i+4, control_allleft_postTFT_dec(i)+8000);
end
%% repeat entire section above for LED trials
led_c1_right_idx=zeros(1,length(tLeftTrial)); %preallocate contrast+side sorted indices
led_c1_left_idx=zeros(1,length(tLeftTrial));
led_c2_right_idx=zeros(1,length(tLeftTrial));
led_c2_left_idx=zeros(1,length(tLeftTrial));
led_c3_right_idx=zeros(1,length(tLeftTrial));
led_c3_left_idx=zeros(1,length(tLeftTrial));
led_c4_right_idx=zeros(1,length(tLeftTrial));
led_c4_left_idx=zeros(1,length(tLeftTrial));

led_c1_righttrials=intersect(find(contrastratio==c1), find(right_led_idx));
led_c1_lefttrials=intersect(find(contrastratio==c1), find(left_led_idx));
led_c1_right_idx(led_c1_righttrials)=1;
led_c1_left_idx(led_c1_lefttrials)=1;
led_c1_rightcorrecttrials=intersect(find(led_c1_right_idx), find(SIx));
led_c1_leftcorrecttrials=intersect(find(led_c1_left_idx), find(SIx));

led_c2_righttrials=intersect(find(contrastratio==c2), find(right_led_idx));
led_c2_lefttrials=intersect(find(contrastratio==c2),find(left_led_idx));
led_c2_right_idx(led_c2_righttrials)=1;
led_c2_left_idx(led_c2_lefttrials)=1;
led_c2_rightcorrecttrials=intersect(find(led_c2_right_idx), find(SIx));
led_c2_leftcorrecttrials=intersect(find(led_c2_left_idx), find(SIx));

led_c3_righttrials=intersect(find(contrastratio==c3), find(right_led_idx));
led_c3_lefttrials=intersect(find(contrastratio==c3), find(left_led_idx));
led_c3_right_idx(led_c3_righttrials)=1;
led_c3_left_idx(led_c3_lefttrials)=1;
led_c3_rightcorrecttrials=intersect(find(led_c3_right_idx), find(SIx));
led_c3_leftcorrecttrials=intersect(find(led_c3_left_idx), find(SIx));

led_c4_righttrials=intersect(find(contrastratio==c4), find(right_led_idx));
led_c4_lefttrials=intersect(find(contrastratio==c4), find(left_led_idx));
led_c4_right_idx(led_c4_righttrials)=1;
led_c4_left_idx(led_c4_lefttrials)=1;
led_c4_rightcorrecttrials=intersect(find(led_c4_right_idx), find(SIx));
led_c4_leftcorrecttrials=intersect(find(led_c4_left_idx), find(SIx));


led_c1_right_plot_qvals=qVals_final(:,led_c1_rightcorrecttrials)-qVals_final(8000,led_c1_rightcorrecttrials);
led_c1_left_plot_qvals=qVals_final(:,led_c1_leftcorrecttrials)-qVals_final(8000,led_c1_leftcorrecttrials);
led_c2_right_plot_qvals=qVals_final(:,led_c2_rightcorrecttrials)-qVals_final(8000,led_c2_rightcorrecttrials);
led_c2_left_plot_qvals=qVals_final(:,led_c2_leftcorrecttrials)-qVals_final(8000,led_c2_leftcorrecttrials);
led_c3_right_plot_qvals=qVals_final(:,led_c3_rightcorrecttrials)-qVals_final(8000,led_c3_rightcorrecttrials);
led_c3_left_plot_qvals=qVals_final(:,led_c3_leftcorrecttrials)-qVals_final(8000,led_c3_leftcorrecttrials);
led_c4_right_plot_qvals=qVals_final(:,led_c4_rightcorrecttrials)-qVals_final(8000,led_c4_rightcorrecttrials);
led_c4_left_plot_qvals=qVals_final(:,led_c4_leftcorrecttrials)-qVals_final(8000,led_c4_leftcorrecttrials);

for i=1:length(led_c1_right_plot_qvals)
led_c1_right_mean_point_qval(i)=nanmean(led_c1_right_plot_qvals(i));
end
for i=1:length(led_c1_left_plot_qvals)
led_c1_left_mean_point_qval(i)=nanmean(led_c1_left_plot_qvals(i));
end
for i=1:length(led_c2_right_plot_qvals)
led_c2_right_mean_point_qval(i)=nanmean(led_c2_right_plot_qvals(i));
end
for i=1:length(led_c2_left_plot_qvals)
led_c2_left_mean_point_qval(i)=nanmean(led_c2_left_plot_qvals(i));
end
for i=1:length(led_c3_right_plot_qvals)
led_c3_right_mean_point_qval(i)=nanmean(led_c3_right_plot_qvals(i));
end
for i=1:length(led_c3_left_plot_qvals)
led_c3_left_mean_point_qval(i)=nanmean(led_c3_left_plot_qvals(i));
end
for i=1:length(led_c4_right_plot_qvals)
led_c4_right_mean_point_qval(i)=nanmean(led_c4_right_plot_qvals(i));
end
for i=1:length(led_c4_left_plot_qvals)
led_c4_left_mean_point_qval(i)=nanmean(led_c4_left_plot_qvals(i));
end
led_mean_qvals_mat=[led_c1_right_mean_point_qval;
    led_c2_right_mean_point_qval;
    led_c3_right_mean_point_qval;
    led_c4_right_mean_point_qval;
    led_c1_left_mean_point_qval;
    led_c2_left_mean_point_qval;
    led_c3_left_mean_point_qval;
    led_c4_left_mean_point_qval;];
colors=brewermap(10, 'RdBu');
colors2=brewermap(10, '*RdBu');

led_c1_rightpostTFT=led_c1_right_mean_point_qval-led_c1_right_mean_point_qval(8800);
led_c1_rightpostTFT_thresh=find(led_c1_rightpostTFT(8800:18000)<=rdt)+8800;
led_c1_rightpostTFT_dec=min(led_c1_rightpostTFT_thresh)-8000;

led_c2_rightpostTFT=led_c2_right_mean_point_qval-led_c2_right_mean_point_qval(8800);
led_c2_rightpostTFT_thresh=find(led_c2_rightpostTFT(8800:18000)<=rdt)+8800;
led_c2_rightpostTFT_dec=min(led_c2_rightpostTFT_thresh)-8000;

led_c3_rightpostTFT=led_c3_right_mean_point_qval-led_c3_right_mean_point_qval(8800);
led_c3_rightpostTFT_thresh=find(led_c3_rightpostTFT(8800:18000)<=rdt)+8800;
led_c3_rightpostTFT_dec=min(led_c3_rightpostTFT_thresh)-8000;

led_c4_rightpostTFT=led_c4_right_mean_point_qval-led_c4_right_mean_point_qval(8800);
led_c4_rightpostTFT_thresh=find(led_c4_rightpostTFT(8800:18000)<=rdt)+8800;
led_c4_rightpostTFT_dec=min(led_c4_rightpostTFT_thresh)-8000; 

led_allright_postTFT_dec=[led_c1_rightpostTFT_dec, led_c2_rightpostTFT_dec, led_c3_rightpostTFT_dec, led_c4_rightpostTFT_dec];
%same for left side below

led_c1_leftpostTFT=led_c1_left_mean_point_qval-led_c1_left_mean_point_qval(8800);
led_c1_leftpostTFT_thresh=find(led_c1_leftpostTFT(8800:18000)>=ldt)+8800;
led_c1_leftpostTFT_dec=min(led_c1_leftpostTFT_thresh)-8000; 

led_c2_leftpostTFT=led_c2_left_mean_point_qval-led_c2_left_mean_point_qval(8800);
led_c2_leftpostTFT_thresh=find(led_c2_leftpostTFT(8800:18000)>=ldt)+8800;
led_c2_leftpostTFT_dec=min(led_c2_leftpostTFT_thresh)-8000;

led_c3_leftpostTFT=led_c3_left_mean_point_qval-led_c3_left_mean_point_qval(8800);
led_c3_leftpostTFT_thresh=find(led_c3_leftpostTFT(8800:18000)>=ldt)+8800;
led_c3_leftpostTFT_dec=min(led_c3_leftpostTFT_thresh)-8000;

led_c4_leftpostTFT=led_c4_left_mean_point_qval-led_c4_left_mean_point_qval(8800);
led_c4_leftpostTFT_thresh=find(led_c4_leftpostTFT(8800:18000)>=ldt)+8800;
led_c4_leftpostTFT_dec=min(led_c4_leftpostTFT_thresh)-8000;

led_allleft_postTFT_dec=[led_c1_leftpostTFT_dec, led_c2_leftpostTFT_dec, led_c3_leftpostTFT_dec, led_c4_leftpostTFT_dec];

for i=1:4
    led_right_actual_dec(i)=led_mean_qvals_mat(i, led_allright_postTFT_dec(i)+8000); %actual quadrature positions at measured RT's
end
for i=1:4
    led_left_actual_dec(i)=led_mean_qvals_mat(i+4, led_allleft_postTFT_dec(i)+8000);
end
colors3=brewermap(10, 'PRGn');
colors4=brewermap(10, '*PRGn');
%
figure 
patch([-2000 0 0 -2000],[-1000 -1000 1000 1000],[0.2 0.2 0.2])
patch([ 0 800 800 0],[-1000 -1000 1000 1000],[0.7 0.7 0.8])
set(gca,'children',flipud(get(gca,'children')))
rdtline=line([-8000 10000], [rdt rdt], 'Color', 'yellow','LineWidth', 3,'LineStyle','--')
ldtline=line([-8000 10000], [ldt ldt], 'Color', 'green', 'LineWidth', 3, 'LineStyle','--')
ylim([-300 400])
xlim([-500 5000])
hold on
p1=plot(qTimes_final,control_c1_right_mean_point_qval,'LineWidth',1.5)
p3=plot(qTimes_final,control_c2_right_mean_point_qval,'LineWidth',1.5 )
p5=plot(qTimes_final,control_c3_right_mean_point_qval,'LineWidth',1.5) 
p7=plot(qTimes_final,control_c4_right_mean_point_qval,'LineWidth',1.5) 
p2=plot(qTimes_final, control_c1_left_mean_point_qval,'LineWidth',1.5 )
p4=plot(qTimes_final, control_c2_left_mean_point_qval,'LineWidth',1.5)
p6=plot(qTimes_final, control_c3_left_mean_point_qval,'LineWidth',1.5)
p8=plot(qTimes_final, control_c4_left_mean_point_qval,'LineWidth',1.5)
p1.Color=colors(1,:)
p2.Color=colors2(1,:)
p3.Color=colors(2,:)
p4.Color=colors2(2,:)
p5.Color=colors(3,:)
p6.Color=colors2(2,:)
p7.Color=colors(4,:)
p8.Color=colors2(4,:)

p9=plot(qTimes_final,led_c1_right_mean_point_qval,'LineWidth',2.5)
p11=plot(qTimes_final,led_c2_right_mean_point_qval,'LineWidth',2.5 )
p13=plot(qTimes_final,led_c3_right_mean_point_qval,'LineWidth',2.5) 
p15=plot(qTimes_final,led_c4_right_mean_point_qval,'LineWidth',2.5) 
p10=plot(qTimes_final, led_c1_left_mean_point_qval,'LineWidth',2.5 )
p12=plot(qTimes_final, led_c2_left_mean_point_qval,'LineWidth',2.5)
p14=plot(qTimes_final, led_c3_left_mean_point_qval,'LineWidth',2.5)
p16=plot(qTimes_final, led_c4_left_mean_point_qval,'LineWidth',2.5)
p9.Color=colors3(1,:)
p10.Color=colors4(1,:)
p11.Color=colors3(2,:)
p12.Color=colors4(2,:)
p13.Color=colors3(3,:)
p14.Color=colors4(2,:)
p15.Color=colors3(4,:)
p16.Color=colors4(4,:)
for i=1:4
        tempvar=find(mean_qvals_mat(i,8000:18000)<= rdt);
         control_right_dectime(i)=min(tempvar);
end
for i=5:8
    tempvar=find(mean_qvals_mat(i,8000:18000)>= ldt);
    control_left_dectime(i-4)=min(tempvar);
end
for i=1:4
        tempvar=find(led_mean_qvals_mat(i,8000:18000)<= rdt);
         led_right_dectime(i)=min(tempvar);
end
for i=5:8
    tempvar=find(led_mean_qvals_mat(i,8000:18000)>= ldt);
    led_left_dectime(i-4)=min(tempvar);
end
inc_resp_check=min(find(mean_qvals_mat(4,8000:10000)>=ldt))+8000; % early wrong answer check
if inc_resp_check<control_c4_rightpostTFT_dec
    control_right_dectime(4)=inc_resp_check
end

rd_yval=[rdt rdt rdt rdt]; ld_yval=[ldt ldt ldt ldt];

led_arerighttimessame=led_allright_postTFT_dec-led_right_dectime;
led_arelefttimessame=led_allleft_postTFT_dec-led_left_dectime;
control_arerighttimessame=control_allright_postTFT_dec-control_right_dectime;
control_arelefttimessame=control_allleft_postTFT_dec-control_left_dectime;
for i=1:4
    if abs(control_arerighttimessame(i))<5
        control_right_dectime(i)=NaN;
    end
    if abs(control_arelefttimessame(i))<5
        control_left_dectime(i)=NaN;
    end
end
for i=1:4
    if abs(led_arerighttimessame(i))<5
        led_right_dectime(i)=NaN;
    end
    if abs(led_arelefttimessame(i))<5
        led_left_dectime(i)=NaN;
    end
end
if input.doSizeDiscrim==1
    taskstr='size';
end
if input.doContrastDiscrim==1
    taskstr='contrast';
end
for i=1:4
    if led_right_dectime(i)>led_allright_postTFT_dec(i) 
        led_right_dectime(i)=NaN;
    end
    if led_left_dectime(i)>led_allleft_postTFT_dec(i)
        led_left_dectime(i)=NaN;
    end
end
for i=1:4
    if control_right_dectime(i)>control_allright_postTFT_dec(i)
        control_right_dectime(i)=NaN;
    end
    if control_left_dectime(i)>control_allleft_postTFT_dec(i)
       control_left_dectime(i)=NaN;
    end
end
if sum(block2)==0
    led_right_dectime=[NaN NaN NaN NaN];
    led_left_dectime=[NaN NaN NaN NaN];
end

scatter(control_right_dectime, rd_yval, 35, 'red','filled','MarkerEdgeColor','black')
scatter(control_left_dectime, ld_yval, 35, 'blue', 'filled','MarkerEdgeColor','black')
scatter(control_allright_postTFT_dec, control_right_actual_dec, 'x','MarkerEdgeColor', 'red', 'LineWidth', 4)
scatter(control_allleft_postTFT_dec, control_left_actual_dec,'x','MarkerEdgeColor','blue','LineWidth',4)
scatter(led_right_dectime, rd_yval, 35, [1 0 1],'filled','MarkerEdgeColor','black')
scatter(led_left_dectime, ld_yval, 35, 'cyan', 'filled','MarkerEdgeColor','black')
scatter(led_allright_postTFT_dec, led_right_actual_dec, 'x','MarkerEdgeColor', [1 0 1], 'LineWidth', 4)
scatter(led_allleft_postTFT_dec, led_left_actual_dec,'x','MarkerEdgeColor','cyan','LineWidth',4)
legend('Stationary period', 'Too fast time', 'Right Decision Threshold', 'Left Decision Threshold', 'Ratio 1 Right', 'Ratio 2 Right Control', 'Ratio 3 Right Control', 'Ratio 4 Right Control','Ratio 1 Left Control', 'Ratio 2 Left Control','Ratio 3 Left Control', 'Ratio 4 Left Control', 'Ratio 1 Right LED','Ratio 2 Right LED', 'Ratio 3 Right LED', 'Ratio 4 Right LED', 'Ratio 1 Left LED','Ratio 2 Left LED','Ratio 3 Left LED', 'Ratio 4 Left LED',  'Right Subjective decision - Control', 'Left Subjective decision - control', 'Right Actual decision - control', 'Left Actual decision - control', 'Right Subjective Decision - LED','Left Subjective Decision - LED','Right Actual Decision - LED', 'Left Actual Decision - LED')   
xlabel('Trial time aligned to stimulus on (ms)'); ylabel('Quadrature Position')
titlestr=sprintf('i435 Average Quadrature Motion Paths by trial side and %s ratio in example session %s', taskstr, savedate(1:6));
title(titlestr);
%% find reaction times with thresholds
uncontrasts=flip(uncontrasts);
figure
hold on
plot(uncontrasts, control_left_dectime','-o','Color', 'blue')
plot(uncontrasts, control_right_dectime, '-o','Color', 'red')
plot(uncontrasts, led_left_dectime,'-o', 'Color', 'cyan')
plot(uncontrasts, led_right_dectime, '-o', 'Color', [1 0 1])
xlabel('Contrast Ratio'); ylabel('Reaction Time (ms)'); title('Reaction Time by contrast ratio and side') 
legend('Left', 'Right')








