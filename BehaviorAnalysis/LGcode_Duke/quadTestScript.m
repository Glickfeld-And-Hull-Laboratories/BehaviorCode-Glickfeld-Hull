Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
cVals_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        cTimes = double(input.counterTimesUs{trN}./1000);
        cVals = double(input.counterValues{trN});
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
                cVals_thresh(:,trN) = cVals(:,find(cTimes>qTimes_thresh(:,trN),1,'first'));
            end
        else
            return
        end
    end
end

minR = input.tooFastTimeMs;
figure;
for it = 1:30
    subplot(5,6,it)
    plot(qTimes_final, qVals_final(:,it))
    xlim([-100 input.tDecisionTimeMs{it}])
    vline(minR)
    if input.tLeftTrial{it}
        title(['Left ' input.trialOutcomeCell{it}])
    else
        title(['Right ' input.trialOutcomeCell{it}])
    end
end


SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'incorrect');
left = cell2mat(input.tLeftTrial);
maxR = input.reactionTimeMs;
lowR = intersect(find(cell2mat(input.tDecisionTimeMs)<maxR), find(cell2mat(input.tDecisionTimeMs)>minR));

qVals_offset = bsxfun(@minus, qVals_final, qVals_final(8001,:));
figure;
subplot(2,2,1)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(SIx),find(left)))))
xlim([-500 2000])
vline(minR)
title('Correct Left Trials')

subplot(2,2,2)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(SIx),find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Correct Right Trials')

subplot(2,2,3)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(MIx),find(left)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Left Trials')

subplot(2,2,4)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(MIx),find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Right Trials')

suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])

print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome.pdf']), '-dpdf')

figure
con_str = strvcat('b', 'r', 'y');
subplot(2,2,1)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(SIx), find(left)))),2), nanstd(qVals_offset(:,(intersect(find(SIx), find(left)))),[],2)./sqrt(length(intersect(find(SIx), find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(SIx), find(left==0)))),2), nanstd(qVals_offset(:,(intersect(find(SIx), find(left==0)))),[],2)./sqrt(length(intersect(find(SIx), find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all correct trials')

subplot(2,2,2)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(MIx), find(left)))),2), nanstd(qVals_offset(:,(intersect(find(MIx), find(left)))),[],2)./sqrt(length(intersect(find(MIx), find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(MIx), find(left==0)))),2), nanstd(qVals_offset(:,(intersect(find(MIx), find(left==0)))),[],2)./sqrt(length(intersect(find(MIx), find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all incorrect trials')

subplot(2,2,3)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(find(SIx), find(left))));
    shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
    hold on;
end
xlim([-500 2000])
vline(minR)
title('Avg left correct trials by contrast')

subplot(2,2,4)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(find(SIx), find(left==0))));
    shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
    hold on
end
xlim([-500 2000])
vline(minR)
title('Avg right correct trials by contrast')
suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome_avg.pdf']), '-dpdf')

data_move = nan(50,nCells,nTrials);
for itrial = 1:nTrials
    if cVals_thresh(itrial)+29 < nframes
        data_move(:,:,itrial) = npSub_tc(cVals_thresh(itrial)-20:cVals_thresh(itrial)+29,:);
    end
end

data_move_dfof = bsxfun(@rdivide, bsxfun(@minus, data_move, dataf), dataf);

figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
       indS = intersect(find(SIx), find(tLeftTrial == iS-1));
       plot(tt,nanmean(data_move_dfof(:,iCell,indS),3));
       hold on;
        if find(good_ind == iCell)
            good_str = ' resp';
        else
            good_str = ' not resp';
        end
        title(['Cell # ' num2str(iCell) ' is' good_str])
    end
end
suptitle([mouse ' ' date '- moveAlign- blue is right; red is left'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_moveAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

