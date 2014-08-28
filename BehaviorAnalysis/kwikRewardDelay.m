function [ output_args ] = kwikRewardDelay(dataFile, startTr, endTr, doBlock)

input = dataFile;
trs = startTr:endTr;

input.holdTimesMs = input.holdTimesMs(trs);
input.reactTimesMs = input.reactTimesMs(trs);
input.juiceTimesMsCell = input.juiceTimesMsCell(trs);
input.tBlock2TrialNumber = input.tBlock2TrialNumber(trs);
input.tDelayTimeMs = input.tDelayTimeMs(trs);
input.holdStartsMs = input.holdStartsMs(trs);
input.trialOutcomeCell = input.trialOutcomeCell(trs);

if doBlock == 1,
    b1Ix = cell2mat_padded(input.tBlock2TrialNumber)==0;
    input.holdTimesMs = input.holdTimesMs(b1Ix);
    input.reactTimesMs = input.reactTimesMs(b1Ix);
    input.juiceTimesMsCell = input.juiceTimesMsCell(b1Ix);
    input.tBlock2TrialNumber = input.tBlock2TrialNumber(b1Ix);
    input.tDelayTimeMs = input.tDelayTimeMs(b1Ix);
    input.holdStartsMs = input.holdStartsMs(b1Ix);
    input.trialOutcomeCell = input.trialOutcomeCell(b1Ix);
end
if doBlock == 2,
    b2Ix = cell2mat_padded(input.tBlock2TrialNumber)==1;
    input.holdTimesMs = input.holdTimesMs(b2Ix);
    input.reactTimesMs = input.reactTimesMs(b2Ix);
    input.juiceTimesMsCell = input.juiceTimesMsCell(b2Ix);
    input.tBlock2TrialNumber = input.tBlock2TrialNumber(b2Ix);
    input.tDelayTimeMs = input.tDelayTimeMs(b2Ix);
    input.holdStartsMs = input.holdStartsMs(b2Ix);
    input.trialOutcomeCell = input.trialOutcomeCell(b2Ix);
end

plotRewardDelay(input, input)



