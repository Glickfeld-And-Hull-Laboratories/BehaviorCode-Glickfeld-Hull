% processing for cerebellarStim(_2P) data that processes errors in trial
% length that point to dropped or added frames and builds an index

dataStruct = struct;
%% formatting assumes "input1" is a cerebellarStim data set
ds = input1;
% and also enforces it
assert(isfield(ds, 'spCounter1TimesUs'), 'Input File is not a cerebellarStim data set');
% takes tr by tr
cTimes = input1.counterTimesUs;
cVals = input1.counterValues;
%  compresses all times/vals
totalTimes = double(cell2mat(cTimes));
totalVals =  double(cell2mat(cVals));
totalTimes = totalTimes(3:end); % removes start-stop-start artifact and 0 value initial state value
totalVals= totalVals(3:end); % removes start-stop-start artifact and 0 value initial state value
%lastTrCounter = cell2mat(input.counter);

imagingRate = input1.frameImagingFrequencyHz; %in Hz
frameIntUs = 1000000/imagingRate;

nTrs = length(input1.counter);
for i=2:nTrs,
    lastVal = cVals{i-1}(end);
    thisTrVal = cVals{i}(1);
    framesBetweenTrs(i-1) = thisTrVal - lastVal;
end
pulseDiff = diff(totalTimes);
pulsePauseIx = pulseDiff>(1.5*frameIntUs);
problematic = pulseDiff(pulsePauseIx);
problemFrameNumbers = find(pulsePauseIx);
normal = pulseDiff(~pulsePauseIx);
avgFrame = mean(normal);

howManyMissed = round((problematic./avgFrame)-1);   %altered by JH 1/14/16 to better reflect how many frames were missed. 
% alert if data has a hole greater than 1 second
assert(sum(howManyMissed>imagingRate)<1, '*** Do not use data: there is a skip of %3d second skip between frames. ****', problematic(howManyMissed>15)/1000000);
hold on
% extract all of our stimulus crap
stimOnFrame = cell2mat(input1.ttCounter);
blank = zeros(1,length(totalVals)+sum(howManyMissed)); %MA 1/20/16: Redefined blank (stimOnIx) to be same length as blankIx, used later. %JH Fix: expands stimOnIx to account for dropped frames. It already adjusts the stimOnFrame to account for missed frames affect on which frame the stimulus corresponds to 
blank(stimOnFrame)=1;
stimOnIx=blank;

%% Check that counter was not reset mid-trial
assert(ismonotonic(totalTimes), 'counterTimesUs is not monotonically increasing: was experiment reset?');
assert(ismonotonic(totalVals), 'counterValues is not monotonically increasing: was experiment reset?');

%% Fill any empty sp/itiCounter10 with the value which was in sp/itiCounter9
for i = 1:input1.stopAfterNTrials
     if isempty(cell2mat(input1.spCounter10(i)))==1
        input1.spCounter10(i) = input1.spCounter9(i);
        sprintf('filled an empty spCounter10')
    end
    if isempty(cell2mat(input1.itiCounter10(i)))==1
        input1.itiCounter10(i) = input1.itiCounter9(i);
        sprintf('filled an empty itiCounter10')
    end
end

%Fill in any missing iti/spCounterTimes
for i = 1:input1.stopAfterNTrials
    if isempty(cell2mat(input1.spCounter10TimesUs(i)))
       input1.spCounter10TimesUs(i) =  mat2cell(cell2mat(input1.spCounter9TimesUs(i))+500000);
        sprintf('filled an empty spCounterTime')
    end
    if isempty(cell2mat(input1.itiCounter10TimesUs(i)))
        input1.itiCounter10TimesUs(i) =  mat2cell(cell2mat(input1.itiCounter9TimesUs(i))+500000);
        sprintf('filled an empty itiCounterTime')
    end
end
    

%% Extracting locomotion data
emptyMat = [];
emptyTimes = [];
locoVals = [];


for i=1:nTrs,
    trN = i;
    
    
    
    itiInt = input1.nScansOff/10;
    spInt = input1.nScansOn/10;
    preInt = input1.preSoundPauseNFrames;
    postInt = input1.postSoundPauseNFrames;
    
    itiMat = ones(1,itiInt);
    spMat = ones(1, spInt);
    stimMat = ones(1,preInt+postInt); 
    
    
    %creates a mat with length = to # of frames per counter and with a
    %value = to the number of wheel pulses in that counter divided by the
    %number of frames in that counter.
    
    iti1Mat = itiMat*double(input1.itiCounter1{trN})/itiInt;
    iti2Mat = itiMat*double(input1.itiCounter2{trN})/itiInt;
    iti3Mat = itiMat*double(input1.itiCounter3{trN})/itiInt;
    iti4Mat = itiMat*double(input1.itiCounter4{trN})/itiInt;
    iti5Mat = itiMat*double(input1.itiCounter5{trN})/itiInt;
    iti6Mat = itiMat*double(input1.itiCounter6{trN})/itiInt;
    iti7Mat = itiMat*double(input1.itiCounter7{trN})/itiInt;
    iti8Mat = itiMat*double(input1.itiCounter8{trN})/itiInt;
    iti9Mat = itiMat*double(input1.itiCounter9{trN})/itiInt;
    iti10Mat = itiMat*double(input1.itiCounter10{trN})/itiInt;
    %Extracts the time each oounter was logged for this trial
    iti1Time = input1.itiCounter1TimesUs{trN}(end);
    iti2Time = input1.itiCounter2TimesUs{trN}(end);
    iti3Time = input1.itiCounter3TimesUs{trN}(end);
    iti4Time = input1.itiCounter4TimesUs{trN}(end);
    iti5Time = input1.itiCounter5TimesUs{trN}(end);
    iti6Time = input1.itiCounter6TimesUs{trN}(end);
    iti7Time = input1.itiCounter7TimesUs{trN}(end);
    iti8Time = input1.itiCounter8TimesUs{trN}(end);
    iti9Time = input1.itiCounter9TimesUs{trN}(end);
    iti10Time = input1.itiCounter10TimesUs{trN}(end);
    %extracts # of wheel pulses in each counter. 
    iti1Val = input1.itiCounter1{trN};
    iti2Val = input1.itiCounter2{trN};
    iti3Val = input1.itiCounter3{trN};
    iti4Val = input1.itiCounter4{trN};
    iti5Val = input1.itiCounter5{trN};
    iti6Val = input1.itiCounter6{trN};
    iti7Val = input1.itiCounter7{trN};
    iti8Val = input1.itiCounter8{trN};
    iti9Val = input1.itiCounter9{trN};
    iti10Val = input1.itiCounter10{trN};

    %repeat process for sp counters
    sp1Mat = spMat*double(input1.spCounter1{trN})/spInt;
    sp2Mat = spMat*double(input1.spCounter2{trN})/spInt;
    sp3Mat = spMat*double(input1.spCounter3{trN})/spInt;
    sp4Mat = spMat*double(input1.spCounter4{trN})/spInt;
    sp5Mat = spMat*double(input1.spCounter5{trN})/spInt;
    sp6Mat = spMat*double(input1.spCounter6{trN})/spInt;
    sp7Mat = spMat*double(input1.spCounter7{trN})/spInt;
    sp8Mat = spMat*double(input1.spCounter8{trN})/spInt;
    sp9Mat = spMat*double(input1.spCounter9{trN})/spInt;
    sp10Mat = spMat*double(input1.spCounter10{trN})/spInt;

    sp1Time = input1.spCounter1TimesUs{trN}(end);
    sp2Time = input1.spCounter2TimesUs{trN}(end);
    sp3Time = input1.spCounter3TimesUs{trN}(end);
    sp4Time = input1.spCounter4TimesUs{trN}(end);
    sp5Time = input1.spCounter5TimesUs{trN}(end);
    sp6Time = input1.spCounter6TimesUs{trN}(end);
    sp7Time = input1.spCounter7TimesUs{trN}(end);
    sp8Time = input1.spCounter8TimesUs{trN}(end);
    sp9Time = input1.spCounter9TimesUs{trN}(end);
    sp10Time = input1.spCounter10TimesUs{trN}(end);

    sp1Val = input1.spCounter1{trN};
    sp2Val = input1.spCounter2{trN};
    sp3Val = input1.spCounter3{trN};
    sp4Val = input1.spCounter4{trN};
    sp5Val = input1.spCounter5{trN};
    sp6Val = input1.spCounter6{trN};
    sp7Val = input1.spCounter7{trN};
    sp8Val = input1.spCounter8{trN};
    sp9Val = input1.spCounter9{trN};
    sp10Val = input1.spCounter10{trN};
    
    
    
    stimMat = stimMat*double(input1.tISIWheelCounter{trN}/(preInt+postInt));   
    stimTime = input1.postSoundWheelCounterTimesUs{trN}(end);
    stimVal = input1.tISIWheelCounter{trN};
    
    %Creates a string with the # of wheel pulses associated with each frame
    %in this trial
    locoMatFrame = [iti1Mat iti2Mat iti3Mat iti4Mat iti5Mat iti6Mat iti7Mat iti8Mat iti9Mat iti10Mat stimMat sp1Mat sp2Mat sp3Mat sp4Mat sp5Mat sp6Mat sp7Mat sp8Mat sp9Mat sp10Mat];
    emptyMat = [emptyMat locoMatFrame]; %begins to concatenate wheel pulse values per frame for all trials 
    locoTimeUs = [iti1Time iti2Time iti3Time iti4Time iti5Time iti6Time iti7Time iti8Time iti9Time iti10Time stimTime sp1Time sp2Time sp3Time sp4Time sp5Time sp6Time sp7Time sp8Time sp9Time sp10Time];
    emptyTimes = [emptyTimes locoTimeUs]; %does the same for counter times. Should not be the same size as locoMatFrames

    tLocoVals =[iti1Val iti2Val iti3Val iti4Val iti5Val iti6Val iti7Val iti8Val iti9Val iti10Val stimVal sp1Val sp2Val sp3Val sp4Val sp5Val sp6Val sp7Val sp8Val sp9Val sp10Val];
    locoVals = [locoVals tLocoVals];
end


%% Replicate graphing to allow frame skipping feedback
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
for i = 1:length(problemFrameNumbers),
    insertFrames = howManyMissed(i);
    startFrame = problemFrameNumbers(i)+1; %to correct for 1 frame "diff" offset
    blankIx(startFrame+shift+1) = 0;     
    emptyMat = cat(2, emptyMat(1:[problemFrameNumbers(i)+shift]), NaN(1,insertFrames), emptyMat([problemFrameNumbers(i)+shift+1]:end));
    stimOnIx((startFrame+shift):(startFrame+insertFrames+shift)) = 0; %make sure this is accurate    DELETE?
    shift = shift+insertFrames;
end
dataStruct.goodFramesIx = boolean(blankIx);
dataStruct.stimOnIx = boolean(stimOnIx)
dataStruct.locomotionMatFrames = emptyMat;
dataStruct.locomotionTimes = emptyTimes;
dataStruct.locomotionVals = locoVals;
%May need to recrop
%%converting locomotionVals and Times to speed in cm/s  ASSUMES
%preSoundPauseNFrames + PostSoundPauseNFrames = num of frames in one iti or
%spCounter

%%WORK IN PROGRESS TO GET SPEED
%frameNumPerBin = diff(dataStruct.locomotionTimes/1000)/250;
%ValuesPerFramedataStruct.locomotionVals(2:end)./frameNumPerBin;
%dataStruct.locomotionSpeed = 

%% CREATING RUN NOTRUN MATRICES AS WELL AS SUSTAINED RUN NOT RUN MATRICES
stillMat = zeros(size(dataStruct.goodFramesIx));
runMat = zeros(size(dataStruct.goodFramesIx));
runThresh = 1.1;
stillThresh =0.9;
for i = 1:size(runMat,2)
    if dataStruct.locomotionMatFrames(i) >= runThresh;
        runMat(i) = 1;
    end
    if dataStruct.locomotionMatFrames(i) < stillThresh;
        stillMat(i) = 1;
    end
end 

 
%% creating matrices for puffs while running and puffs while still

puffRun = and(runMat,dataStruct.stimOnIx);
puffStill = and(stillMat, dataStruct.stimOnIx);
dataStruct.puffRun = puffRun;
dataStruct.puffStill = puffStill;


%% WRITE CODE HERE TO CROP OUT TRIALS WITH 1 SECOND HOLES

%% WRITE CODE HERE TO IMPORT IMAGE
%dataStruct.img = imread('');
%% WRITE CODE HERE TO INDEX IMAGE

%% WRITE CODE HERE TO CHECK IF IMAGE HAS THE SAME NUMBER OF FRAMES AS BEHAVIOR FILE HAS IN COUNTERS

%% IF doLED ACTIVE, SEARCH FOR HIGH LUMINANCE SPIKES

% if isfield(input1, 'doLED')
%     if input1.doLED==1,
    % do LED searchy type things: probably something like "all the luminous frames please stand up"
    % and then correlate that with trial starts/ends to see if you're all
    % matched up
%     end
% end
