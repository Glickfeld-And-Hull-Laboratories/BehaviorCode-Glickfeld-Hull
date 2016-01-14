power = cell2mat(input.tTrialLaserPowerMw);
outcomes = input.trialOutcomeCell;
nLevels = unique(power);
subplotSz = {2,1};
for x=1:length(nLevels);
    i=nLevels(x);
    iX = power==i;
    levelOutcomes = outcomes(iX);
    nTrs(x) = length(levelOutcomes);
    successIx = strcmp(levelOutcomes, 'success');
    perCorrect(x) = sum(successIx)/nTrs(x);
    earlyIx = strcmp(levelOutcomes, 'failure');
    perEarly(x) = sum(earlyIx)/nTrs(x);
    perIgnore(x) = sum(strcmp(levelOutcomes, 'ignore'))/nTrs(x);
    
    holds = cell2mat(input.holdTimesMs(iX));
    meanCHolds(x) = mean(holds(successIx));
    meanEHolds(x) = mean(holds(earlyIx));
    
    
end

%subplot(subplotSz{:},1);
hold on
plot(nLevels, perCorrect, 'o-k')
plot(nLevels, perEarly, 'o-c')
plot(nLevels, perIgnore, 'o-m')
ylim([0 1])
grid on
legend('Correct Trials', 'Early Trials', 'Ignored Trials')

%subplot(subplotSz{:},2);
%hold on
%plot(nLevels, meanCHolds, 'o-k');
%plot(nLevels, meanEHolds, 'o-c');
