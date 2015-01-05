% So I've added a little script to flashingStimTest.m that will create a new outcome field: catchTrialOutcomeCell.
% 'NaN'- means it was not a short catch trial
% 'FA' - means it was a false alarm (the mouse released in response to the catch stimulus)
% 'CR' - means it was a correct reject (the mouse held through the react period)
% 'failure' - means that the mouse did not hold long enough to reach the catch stimulus

for trN = 1:length(missIx)
  if input.tShortCatchTrial{trN}
    if input.tFalseAlarm{trN}
        input.catchTrialOutcomeCell{trN} = 'FA';
    end
    if isempty(input.tCatchTimeMs{trN})
        input.tCatchTimeMs{trN} = NaN;
        input.catchTrialOutcomeCell{trN} = 'failure';
    end
    if (input.leverUpTimeMs{trN}-input.tCatchTimeMs{trN})>input.reactTimeMs
        input.catchTrialOutcomeCell{trN} = 'CR';
    end
    if (input.leverUpTimeMs{trN}-input.tCatchTimeMs{trN})<input.tooFastTimeMs
        input.catchTrialOutcomeCell{trN} = 'failure';
    end
else
    input.catchTrialOutcomeCell{trN} = 'NaN';
end;
end