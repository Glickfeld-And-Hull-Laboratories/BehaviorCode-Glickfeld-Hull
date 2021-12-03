function tLeftChoice = calcChoice2AFC(input);
%uses quadrature output to recalculate direction of choice
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


