function tLeftChoice = calcChoice2AFC(input);
%uses quadrature output to recalculate direction of choice
if isfield(input,'counterValues')
    qReact = celleqel2mat_padded(input.qStartReact);
    cDecision = celleqel2mat_padded(input.cDecision);
    nTrials = size(cDecision,2);
    qDecision = zeros(size(cDecision));
    for i = 1:nTrials
         counterVals = input.counterValues{i};
         counterTimes = input.counterTimesUs{i};
         quadVals = input.quadratureValues{i};
         quadTimes = input.quadratureTimesUs{i};
         tDecision = counterTimes(find(counterVals == cDecision(i)));
         [min_val min_ind] = min(abs(quadTimes-tDecision));
         qDecision(i) = quadVals(min_ind);
    end
    tLeftChoice = qDecision-qReact > 0;
else
    qReact = celleqel2mat_padded(input.qStartReact);
    decisionTime = celleqel2mat_padded(input.tDecisionTimeMs);
    nTrials = size(decisionTime,2);
    qDecision = zeros(size(decisionTime));
    for i = 1:nTrials
         quadVals = input.quadratureValues{i};
         quadTimes = input.quadratureTimesUs{i};
         if length(input.delayTimeMs)>1
             input.delayTimeMs = input.delayTimeMs(1);
         end
         tDecision = input.stimTimestampMs{i}+decisionTime(i)-input.delayTimeMs;
         [min_val min_ind] = min(abs(quadTimes-(tDecision*1000)));
         qDecision(i) = quadVals(min_ind);
    end
    tLeftChoice = qDecision-qReact > 0;
end

