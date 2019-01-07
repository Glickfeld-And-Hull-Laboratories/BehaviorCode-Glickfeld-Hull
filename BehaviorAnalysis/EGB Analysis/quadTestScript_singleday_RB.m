cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
%cVals_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
stimtime_vector=input.stimTimestampMs;
tContrast=celleqel2mat_padded(input.tGratingContrast);
dContrast=celleqel2mat_padded(input.dGratingContrast);
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
%         cTimes = double(input.counterTimesUs{trN}./1000);
%         cVals = double(input.counterValues{trN});
        stimTime = double(input.stimTimestampMs{trN});
        %stimVal = double(input.qStimOn{trN});
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
               % cVals_thresh(:,trN) = cVals(:,find(cTimes>qTimes_thresh(:,trN),1,'first'));
            end
        else
            return
        end
    end
end


%%
ldt=input.leftDecisionThreshold;
rdt=input.rightDecisionThreshold;
c1_righttrials=intersect(find(contrastratio==c1), find(~tLeftTrial));
c1_lefttrials=intersect(find(contrastratio==c1), find(tLeftTrial));
c2_righttrials=intersect(find(contrastratio==c2), find(~tLeftTrial));
c2_lefttrials=intersect(find(contrastratio==c2),find(tLeftTrial));
c3_righttrials=intersect(find(contrastratio==c3), find(~tLeftTrial));
c3_lefttrials=intersect(find(contrastratio==c3),find(tLeftTrial));
c4_righttrials=intersect(find(contrastratio==c4), find(~tLeftTrial));
c4_lefttrials=intersect(find(contrastratio==c4),find(tLeftTrial));

c1_right_plot_qtimes= qTimes_act(:,c1_righttrials)-stimtime_vector(c1_righttrials);
c1_left_plot_qtimes=qTimes_act(:,c1_lefttrials)-stimtime_vector(c1_lefttrials);

c1_right_plot_qvals=qVals_final(:,c1_righttrials)-qVals_final(8000,c1_righttrials);
c1_left_plot_qvals=qVals_final(:,c1_lefttrials)-qVals_final(8000,c1_lefttrials);
c2_right_plot_qvals=qVals_final(:,c2_righttrials)-qVals_final(8000,c2_righttrials);
c2_left_plot_qvals=qVals_final(:,c2_lefttrials)-qVals_final(8000,c2_lefttrials);
c3_right_plot_qvals=qVals_final(:,c3_righttrials)-qVals_final(8000,c3_righttrials);
c3_left_plot_qvals=qVals_final(:,c3_lefttrials)-qVals_final(8000,c3_lefttrials);
c4_right_plot_qvals=qVals_final(:,c4_righttrials)-qVals_final(8000,c4_righttrials);
c4_left_plot_qvals=qVals_final(:,c4_lefttrials)-qVals_final(8000,c4_lefttrials);

for i=1:length(c1_right_plot_qvals)
c1_right_mean_point_qval(i)=nanmean(c1_right_plot_qvals(i));
end
for i=1:length(c1_left_plot_qvals)
c1_left_mean_point_qval(i)=nanmean(c1_left_plot_qvals(i));
end
for i=1:length(c2_right_plot_qvals)
c2_right_mean_point_qval(i)=nanmean(c2_right_plot_qvals(i));
end
for i=1:length(c2_left_plot_qvals)
c2_left_mean_point_qval(i)=nanmean(c2_left_plot_qvals(i));
end
for i=1:length(c3_right_plot_qvals)
c3_right_mean_point_qval(i)=nanmean(c3_right_plot_qvals(i));
end
for i=1:length(c3_left_plot_qvals)
c3_left_mean_point_qval(i)=nanmean(c3_left_plot_qvals(i));
end
for i=1:length(c4_right_plot_qvals)
c4_right_mean_point_qval(i)=nanmean(c4_right_plot_qvals(i));
end
for i=1:length(c4_left_plot_qvals)
c4_left_mean_point_qval(i)=nanmean(c4_left_plot_qvals(i));
end

figure 
hold on
%plot(c1_right_plot_qtimes, c1_right_plot_qvals,'Color',[ 0 1 1]) %cyan
%plot(c1_left_plot_qtimes, c1_left_plot_qvals,'Color',[1 0.5 0]) %orange
patch([-2000 0 0 -2000],[-200 -200 400 400],[0.7 0.7 0.8])
patch([ 0 800 800 0],[-200 -200 400 400],[0.2 0.2 0.2])
set(gca,'children',flipud(get(gca,'children')))

ylim([-200 400])
xlim([-2000 10001])
plot(qTimes_final,c1_right_mean_point_qval, 'Color', 'blue','LineWidth', 1.5, 'LineStyle', '--')
plot(qTimes_final, c1_left_mean_point_qval, 'Color', 'red', 'Linewidth', 1.5)
plot(qTimes_final,c2_right_mean_point_qval, 'Color', [0 0 0.85],'LineWidth',1.5, 'LineStyle', '--')
plot(qTimes_final, c2_left_mean_point_qval, 'Color', 'cyan', 'Linewidth', 1.5)
plot(qTimes_final,c3_right_mean_point_qval, 'Color', [0 0 0.7],'LineWidth', 1.5, 'LineStyle', '--')
plot(qTimes_final, c3_left_mean_point_qval, 'Color', [0.5 1 0.2], 'Linewidth', 1.5)
plot(qTimes_final,c4_right_mean_point_qval, 'Color', [0 0 0.5],'LineWidth', 1.5, 'LineStyle', '--')
plot(qTimes_final, c4_left_mean_point_qval, 'Color', [0 1 1], 'Linewidth', 1.5)
rdtline=line([-8000 10000], [-rdt -rdt], 'Color', 'yellow','LineWidth', 3,'LineStyle','--')
ldtline=line([-8000 10000], [ldt ldt], 'Color', 'yellow', 'LineWidth', 3, 'LineStyle','--')
%legend(' ', '', '', 'AVG right trial motion', 'AVG left trial motion', 'All right trials', 'All left trials')
title('i414 quadrature motion paths in example session'); xlabel('Time (ms)'); ylabel('Quadrature position')
%% find reaction times with thresholds








