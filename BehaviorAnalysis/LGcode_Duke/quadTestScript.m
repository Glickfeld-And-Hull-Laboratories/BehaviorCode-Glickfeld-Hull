Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        stimTime = double(input.stimTimestampMs{trN});
        %stimVal = double(input.qStimOn{trN});
        qTimes = qTimes-stimTime;
        qVals = qVals-qVals(1);
        time_ind = find(qTimes>= -8000 & qTimes<=10000);
        if length(time_ind>2)
            qTimes_sub = qTimes(time_ind);
            qVals_sub = qVals(time_ind);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qTimes_final = -8000:10000;
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)';
        else
            return
        end
    end
end

figure;
for it = 1:length(ind)
    trN = x(ind(it));
    subplot(5,6,it)
    plot(qTimes_final, qVals_final(:,trN))
    xlim([-100 input.tDecisionTimeMs{trN}])
    vline(750)
    if input.tLeftTrial{trN}
        title(['Left ' input.trialOutcomeCell{trN}])
    else
        title(['Right ' input.trialOutcomeCell{trN}])
    end
end


SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'incorrect');
left = cell2mat(input.tLeftTrial);
minR = 750;
maxR = 25000;
lowR = intersect(find(cell2mat(input.tDecisionTimeMs)<maxR), find(cell2mat(input.tDecisionTimeMs)>minR));

figure;
subplot(2,2,1)
plot(qTimes_final, bsxfun(@minus, qVals_final(:,intersect(lowR, intersect(find(SIx),find(left)))),qVals_final(8001,intersect(lowR, intersect(find(SIx),find(left))))))
xlim([-500 2000])
vline(750)
title('Correct Left Trials')

subplot(2,2,2)
plot(qTimes_final, bsxfun(@minus, qVals_final(:,intersect(lowR, intersect(find(SIx),find(left==0)))),qVals_final(8001,intersect(lowR, intersect(find(SIx),find(left==0))))))
xlim([-500 2000])
vline(750)
title('Correct Right Trials')

subplot(2,2,3)
plot(qTimes_final, bsxfun(@minus, qVals_final(:,intersect(lowR, intersect(find(MIx),find(left)))),qVals_final(8001,intersect(lowR, intersect(find(MIx),find(left))))))
xlim([-500 2000])
vline(750)
title('Incorrect Left Trials')

subplot(2,2,4)
plot(qTimes_final, bsxfun(@minus, qVals_final(:,intersect(lowR, intersect(find(MIx),find(left==0)))),qVals_final(8001,intersect(lowR, intersect(find(MIx),find(left==0))))))
xlim([-500 2000])
vline(750)
title('Incorrect Right Trials')

suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])

print(['Z:\home\lindsey\Analysis\Behavior\2AFC_Quadrature\' input.saveTime '-' num2str(input.subjectNum) '-ReactRange-'  num2str(minR) '-' num2str(maxR) 'ms.pdf'],'-dpdf')

figure
subplot(2,2,1)
plot(qTimes_final, nanmean(bsxfun(@minus,qVals_final(:,intersect(lowR,intersect(find(SIx), find(left)))), qVals_final(8001,intersect(lowR,intersect(find(SIx), find(left))))),2))
hold on;
plot(qTimes_final, nanmean(bsxfun(@minus,qVals_final(:,intersect(lowR,intersect(find(SIx), find(left==0)))), qVals_final(8001,intersect(lowR,intersect(find(SIx), find(left==0))))),2))
xlim([-500 2000])
vline(750)
title('Avg all correct trials')

subplot(2,2,2)
plot(qTimes_final, nanmean(bsxfun(@minus,qVals_final(:,(intersect(find(MIx), find(left)))), qVals_final(8001,(intersect(find(MIx), find(left))))),2))
hold on;
plot(qTimes_final, nanmean(bsxfun(@minus,qVals_final(:,(intersect(find(MIx), find(left==0)))), qVals_final(8001,(intersect(find(MIx), find(left==0))))),2))
xlim([-500 2000])
vline(750)
title('Avg all incorrect trials')
print(['Z:\home\lindsey\Analysis\Behavior\2AFC_Quadrature\' input.saveTime '-' num2str(input.subjectNum) '-ReactRange-'  num2str(minR) '-' num2str(maxR) 'ms-AVG.pdf'],'-dpdf')
