subNum = '613';
date = '150508';
time = '1430';

CD = 'Y:\home\andrew\Behavior\Data';
cd(CD);

mworks = ['data-i' subNum '-' date '-' time '.mat'];
load(mworks)

%%
% datasets without input.catchTrialOutcomeCell
nTrials = input.trialSinceReset;


for trN = 1:nTrials
  if input.tShortCatchTrial{trN}
      if strcmp(input.trialOutcomeCell{trN},'success') == 1
          input.catchTrialOutcomeCell{trN} = 'CR';
      else 
          input.catchTrialOutcomeCell{trN} = 'failure';
      end
      if input.tFalseAlarm{trN}
          input.catchTrialOutcomeCell{trN} = 'FA';
      end
  else
        input.catchTrialOutcomeCell{trN} = 'NaN';
  end;
end

%%

DirectionDeg = cell2mat_padded(input.tGratingDirectionDeg);
doCatch = cell2mat_padded(input.tCatchGratingDirectionDeg);
catchDeg = unique(doCatch);
trialOutcome = input.trialOutcomeCell;
catchOutcome = input.catchTrialOutcomeCell;
doBlock2 = cell2mat_padded(input.tBlock2TrialNumber);
visTrialOutcome = trialOutcome(find(doBlock2 == 0));

% %last half of training
% nTrials = input.trialSinceReset;
% midTrial = nTrials/2;
% DirectionDeg = DirectionDeg(midTrial:end);
% doCatch = doCatch(midTrial:end);
% trialOutcome = trialOutcome(midTrial:end);
% catchOutcome = catchOutcome(midTrial:end);
% doBlock2 = doBlock2(midTrial:end);
% visTrialOutcome = trialOutcome(find(doBlock2 == 0));

% deg of change for two levels of catch trial
catch2 = catchDeg(end);
catch1 = catchDeg(2);

% find catch trials for specific level of change, this is the total catch
% trials for each level
catchTrials1 = find(doCatch == catch1);
catchTrials2= find(doCatch == catch2);

%find outcome for those catch trials, FA means mouse released in response
%to the catch stimulus (visual change).
outcomeCatch1 = catchOutcome(catchTrials1);
outcomeCatch2 = catchOutcome(catchTrials2);

%find visual trials for same degrees of change as catch stimulus
visTrials1 = find(DirectionDeg == catch1 & doBlock2 == 0);
visTrials2 = find(DirectionDeg == catch2 & doBlock2 == 0);

%find outcome for those visual trials; success means mouse release in
%response to visual change.
outcomeVis1 = trialOutcome(visTrials1);
outcomeVis2 = trialOutcome(visTrials2);

%find percentage of responses to catch stim and vis stim
respYesCatchTrials1 = length(find(strcmp('FA',outcomeCatch1)));
respYesCatchTrials2 = length(find(strcmp('FA',outcomeCatch2)));

percentYesCatchTrials1 = (respYesCatchTrials1/length(outcomeCatch1))*100;
percentYesCatchTrials2 = (respYesCatchTrials2/length(outcomeCatch2))*100;
percentYesCatchMat = [percentYesCatchTrials1 percentYesCatchTrials2]';

respYesVisTrials1 = length(find(strcmp('success',outcomeVis1)));
respYesVisTrials2 = length(find(strcmp('success',outcomeVis2)));

percentYesVisTrials1 = (respYesVisTrials1/length(outcomeVis1))*100;
percentYesVisTrials2 = (respYesVisTrials2/length(outcomeVis2))*100;
percentYesVisMat = [percentYesVisTrials1 percentYesVisTrials2]';

%bar graph
catchExpResults = cat(2,percentYesCatchMat,percentYesVisMat);
legendinfo = {'catch trials', ; 'visual trials'}
figure;
bar(catchExpResults,'grouped')
hold on
legend(legendinfo)
title(['short catch trials - ' date])
ylabel('percent respond yes')
set(gca,'XTick',1:2,'XTickLabel',[catch1 catch2])

