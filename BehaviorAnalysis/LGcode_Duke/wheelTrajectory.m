function [qVals_final qTimes_thresh] = wheelTrajectory(input,thresh);
%qVals_final gives wheel position from -8s to 10s from stim on
%qTimes_thresh gives the time that the wheel trajectory is >thresh ticks

Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));

if isfield(input, 'didNoGo')
    Nix = find(celleqel2mat_padded(input.didNoGo));
    Tix = setdiff(1:length(input.trialOutcomeCell), [Nix Iix]);
else
    Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
end
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
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
                qTimes_thresh(:,trN) = find(abs(qVal_off(8000:end,:))>thresh,1,'first');
            end
        else
            return
        end
    end
end