%New cerebellarFrameAnalyzer for updated data collection program.

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
totalTimes = totalTimes(3:end); %removes start-stop-start artifact.
totalVals = totalVals(3:end); %removes start-stop-start artifact.

frameIntMs =  double(input1.frameImagingDurationMs);
imagingRate = 1000/frameIntMs;

pulseDiff = diff(totalTimes);
pulsePauseIx = pulseDiff>(1.5*frameIntMs);

%create index for when stimulus was delivered
stimOnFrame = cell2mat(input1.tStimTurnedOn);
stimOnIx = zeros(1, length(totalTimes));
for ii=1:length(stimOnFrame)
    stimOn = find((totalTimes/stimOnFrame(ii)) < 1.00001 & (totalTimes/stimOnFrame(ii))> 0.99999);  %creates a logical array that indexes which frames were stimulus frames, according to their timestamp.
    stimOnIx(stimOn) = 1;
end

%% Check that counter was not reset mid-trial
assert(ismonotonic(totalTimes), 'Counter times are not monotonically increasing. Was experiment reset?');
assert(ismonotonic(totalVals), ' Counter values are not monotonically increasing. Was experiment reset?');
%% Extract and format locomotion data

quadTimes = double(cell2mat(input1.quadratureTimesUs))/1000;
quadValues = double(cell2mat(input1.quadratureValues));
%Removes start-stop-start artifact data from quadTimes/quadValues
startquad = find((quadTimes/totalTimes(1)) < 1.000001 & (quadTimes/totalTimes(1)) > 0.999999);
quadTimes = quadTimes(startquad:end); quadValues = quadValues(startquad:end);

%Generates an index in which each value points to the location of the start
%of a frame in quadTimes and quadValues.
quadFrameIndex = NaN(1,length(totalTimes));
for ii=1:length(totalTimes)
    quadFrame = find((quadTimes/totalTimes(ii)) < 1.000001 & (quadTimes/totalTimes(ii)) > 0.999999);
    quadFrameIndex(:,ii) = quadFrame;
end

%Creates a frame-by-frame index of wheel pulse data (pulses per frame).
locoMat = NaN(1,length(totalTimes)-1);
for ii=2:(length(quadFrameIndex))
    locoMatFrame = sum(diff(quadValues(quadFrameIndex(ii-1):quadFrameIndex(ii))));
    locoMat(:,ii-1) = locoMatFrame;
end
locoMat = [locoMat locoMat(end)];

%% Replicate graphing to allow frame skipping feedback
%Copied from original cerebellarFrameAnalyzer to check for dropped frames.
problematic = pulseDiff(pulsePauseIx);
problemFrameNumbers = find(pulsePauseIx);
normal = pulseDiff(~pulsePauseIx);
avgFrame = mean(normal);
howManyMissed = (problematic./avgFrame)-1;

if length(problemFrameNumbers)>0,
    [a b c] = plotyy(1:length(pulseDiff), pulseDiff, problemFrameNumbers, howManyMissed)
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
    locoMat = cat(2, locoMat(1:[problemFrameNumbers(i)+shift]), NaN(1,insertFrames), locoMat([problemFrameNumbers(i)+shift+1]:end));
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
