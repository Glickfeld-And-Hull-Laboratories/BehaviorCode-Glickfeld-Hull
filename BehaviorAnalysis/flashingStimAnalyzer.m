%New cerebellarFrameAnalyzer for visual flashing stimulus.

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

%create index for when stimulus was delivered
input1.movieStartMS{1} = double(input1.movieStartMS{1});
flashTimes = zeros(1,200);
flashTimes(1) = input1.movieStartMS{1}+11000;       %Movie starts before first frame, so the first "real" stimulus is the second one in the movie. (11 seconds between the start of each flash)
for ii=2:200
    flashTimes(ii) = flashTimes(ii-1)+11000;
end
stimOnIx = zeros(1, length(totalTimes)+sum(howManyMissed));
stimStartStopIx = zeros(1,length(totalTimes) + sum(howManyMissed));
for ii=1:length(flashTimes)
    flashOn = find(totalTimes<flashTimes(ii),1,'last');  %Finds the frame at which the visual movie began in order to index all subsequent stimuli
    if flashTimes(ii)>totalTimes(end)
        break
    end
    stimOnIx(flashOn) = 1;
    stimStartStopIx(flashOn) = 1; stimStartStopIx(flashOn+10)=1;
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
runThresh = 5;
stillThresh = 4;
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