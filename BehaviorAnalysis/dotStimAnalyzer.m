%New cerebellarFrameAnalyzer for moving dots stimulus.

dataStruct = struct;
%% formatting assumes "input1" is a cerebellarStim data set
ds = input1;
assert(isfield(ds,'quadratureTimesUs'), 'Input File is not a CerebellarStim data set');
%takes trial by trial
cTimes = input1.counterTimesUs;
cVals = input1.counterValues;
%compresses all times/vals
totalTimes = double(cell2mat(cTimes))/1000; %Converts counterTimes to ms.
totalVals = double(cell2mat(cVals));
%removes start-stop-start artifact (assumes last value of counterValues is
%largest value, which is almost always the case)
startframe = (length(totalVals)-max(totalVals))+1;
totalTimes = totalTimes(startframe:end);
totalVals = totalVals(startframe:end); 

frameIntMs =  double(input1.frameImagingDurationMs);
imagingRate = 1000/frameIntMs;

pulseDiff = diff(totalTimes);
pulsePauseIx = pulseDiff>(1.5*frameIntMs);

problematic = pulseDiff(pulsePauseIx);
problemFrameNumbers = find(pulsePauseIx);
normal = pulseDiff(~pulsePauseIx);
avgFrame = mean(normal);
howManyMissed = round((problematic./avgFrame)-1);

%create index for when stimulus was delivered (brake)
if input1.DoSolenoid==1
    brake_times = input1.tResistanceStartMS; brake_times = cell2mat(brake_times);
    if brake_times(1)==0;
        brake_times(1)=[];
    end
    brake_dur = input1.solenoidStimulusDurationMs/100;
    stimOnIx = zeros(1,length(totalTimes)+sum(howManyMissed));
    stimStartStopIx = zeros(1,length(totalTimes)+sum(howManyMissed));
    for ii=1:length(brake_times)
        brake_on = find(totalTimes<brake_times(ii),1,'last');
        if brake_times(ii)>totalTimes(end)
            break
        end
        stimOnIx(brake_on:brake_on+brake_dur-1) = 1;
        stimStartStopIx(brake_on) = 1; stimStartStopIx(brake_on+brake_dur-1) = 1;
    end
end    

%create index for when stimulus was delivered (optical flow reversal)
if input1.DoMovingDot==1 && input1.DoCerebellarStim==1
    reverse_times = input1.tReverseVStimTimeMs; reverse_times = cell2mat(reverse_times);
    if reverse_times(1)==0;
        reverse_times(1)=[];
    end
    reverse_dur = input1.ReverseDurationMS/100;
    stimOnIx = zeros(1,length(totalTimes)+sum(howManyMissed));
    stimStartStopIx = zeros(1,length(totalTimes)+sum(howManyMissed));
    for ii=1:length(reverse_times)
        rev_on = find(totalTimes<reverse_times(ii),1,'last');
        if reverse_times(ii)>totalTimes(end)
            break
        end
        stimOnIx(rev_on:rev_on+reverse_dur-1) = 1;
        stimStartStopIx(rev_on) = 1; stimStartStopIx(rev_on+reverse_dur-1) = 1;
    end
end

%create index for when stimulus was delivered (moving dots)
if input1.DoMovingDot==1 && input1.DoCerebellarStim==0;
    visOn = cell2mat(input1.tVisualStimTurnedOnMs);
    visOff = cell2mat(input1.tVisualStimTurnedOffMs);
    visDuration = (visOff./1000) - (visOn./1000); visDuration = round(visDuration.*10);

    stimOnIx = zeros(1, length(totalTimes)+sum(howManyMissed));
    stimStartStopIx = zeros(1, length(totalTimes) + sum(howManyMissed));
    for ii=1:length(visOn)
        dotsOn = find(totalTimes<visOn(ii),1,'last');  
        if visOn(ii)>totalTimes(end)
            break
        end
        stimOnIx(dotsOn:dotsOn+visDuration(ii)-1) = 1;
        stimStartStopIx(dotsOn) = 1; stimStartStopIx(dotsOn+visDuration(ii)-1)=1;
    end
end

%% Check that counter was not reset mid-trial
assert(ismonotonic(totalTimes), 'Counter times are not monotonically increasing. Was experiment reset?');
assert(ismonotonic(totalVals), ' Counter values are not monotonically increasing. Was experiment reset?');
%% Extract and format locomotion data
%Error fix: Empty cells in quadValues/quadTimes (mouse didn't move for that
%trial) have different data types (double) than the other cells (int64).
%This loop finds those cells and converts them to int64 so we can use
%cell2mat below. MA 2/3/16
if sum(cellfun('isempty',input1.quadratureValues))> 0
    emptyCellIx = find(cellfun('isempty',input1.quadratureValues));
    for ii=1:length(emptyCellIx)
        input1.quadratureValues{emptyCellIx(ii)}=int64(input1.quadratureValues{emptyCellIx(ii)});
        input1.quadratureTimesUs{emptyCellIx(ii)}=int64(input1.quadratureTimesUs{emptyCellIx(ii)});
    end
end

quadTimes = double(cell2mat(input1.quadratureTimesUs))/1000;
quadValues = double(cell2mat(input1.quadratureValues));
%Removes start-stop-start artifact data from quadTimes/quadValues
startquad = find(quadTimes>totalTimes(1),1,'first');
if isempty(startquad)==1
    startquad = find(quadTimes<totalTimes(1),1,'first');
    quadTimes = quadTimes(startquad:end); quadValues = quadValues(startquad:end);
end
quadTimes = quadTimes(startquad:end); quadValues = quadValues(startquad:end);
quadDiff = diff(quadValues); quadDiff = [1 quadDiff];
quadInd = zeros(1, length(quadTimes));
locoMat = zeros(1, length(totalTimes));

for ii=1:length(quadTimes)
    quadInd(ii) = find(totalTimes<quadTimes(ii),1,'last');
    locoMat(quadInd(ii)) = locoMat(quadInd(ii))+ quadDiff(ii);
end

%% Replicate graphing to allow frame skipping feedback
%Copied from original cerebellarFrameAnalyzer to check for dropped frames.
if isempty(problemFrameNumbers)==0,
    [a, b, c] = plotyy(1:length(pulseDiff), pulseDiff, problemFrameNumbers, howManyMissed);
    title('Estimates of Missed Frames')
    set(a(1), 'Ylim', [0 max(pulseDiff)])
    set(a(2), 'Ylim', [0 max(howManyMissed)])
    set(c, 'Marker', 'x')
    set(c, 'LineStyle', 'none')
    ylabel('Inter-Frame Latency (us)')
    xlabel('Frames')
    set(get(a(2),'Ylabel'),'String','Frames Missed')
else
    plot(1:length(pulseDiff), pulseDiff)
    ylim([0 max(pulseDiff)])
    xlabel('Frames')
    ylabel('Inter-Frame Latency (us)')
end

%% Generate the good frames index
shift = 0;
blankIx = ones(1,length(totalVals)+sum(howManyMissed));
for i = 1:length(problemFrameNumbers)
    insertFrames = howManyMissed(i);
    startFrame = problemFrameNumbers(i)+1; %to correct for 1 frame "diff" offset
    blankIx(startFrame+shift+1) = 0;     
    locoMat = cat(2, locoMat(1:problemFrameNumbers(i)+shift), NaN(1,insertFrames), locoMat(problemFrameNumbers(i)+shift+1:end));
    shift = shift+insertFrames;
end
dataStruct.goodFramesIx = boolean(blankIx);
dataStruct.stimOnIx = boolean(stimOnIx);
dataStruct.locomotionMatFrames = locoMat;

%% Create run/not-run matrices
%Copied from old cerebellarFrameAnalyzer code, except now compare absolute
%value of locoMat values due to possibility of negative velocity.
stillMat = zeros(size(dataStruct.goodFramesIx));
runMat = zeros(size(dataStruct.goodFramesIx));
runThresh = 3;
stillThresh = 2;
for ii = 1:length(runMat)
    if abs(dataStruct.locomotionMatFrames(ii)) >= runThresh;
        runMat(ii) = 1;
    end
    if abs(dataStruct.locomotionMatFrames(ii)) < stillThresh;
        stillMat(ii) = 1;
    end
end
puffRun = and(runMat,dataStruct.stimOnIx);
puffStill = and(stillMat,dataStruct.stimOnIx);
dataStruct.puffRun = puffRun;
dataStruct.puffStill = puffStill;
dataStruct.stimStartStopIx = stimStartStopIx;
