% function [dataStruct] = puffNRun(folder)
pathName = 'Z:\Data\WidefieldImaging\GCaMP\150409_img22_1';
imgName = [pathName '\150409_img22_1_MMStack.ome'];
input = load([pathName '\data-img22-150409-1621']);
if isfield(input,'input'),
    input = input.input;
end
info = imfinfo(imgName, 'TIF');


%%Fill any empty sp/itiCounters with a NaN
for i = 1:input.stopAfterNTrials
     if isempty(cell2mat(input.spCounter10(i)))==1
        input.spCounter10(i) = input.spCounter9(i);
        sprintf('filled an empty spCounter10')
    end
    if isempty(cell2mat(input.itiCounter10(i)))==1
        input.itiCounter10(i) = input.itiCounter9(i);
        sprintf('filled an empty itiCounter10')
    end
end

%Fill in any missing iti/spCounterTimes
for i = 1:input.stopAfterNTrials
    if isempty(cell2mat(input.spCounter10TimesUs(i)))
        input.spCounter10TimesUs(i) =  mat2cell(cell2mat(input.spCounter9TimesUs(i))+500000);
        sprintf('filled an empty spCounterTime')
    end
    if isempty(cell2mat(input.itiCounter10TimesUs(i)))
        input.itiCounter10TimesUs(i) =  mat2cell(cell2mat(input.itiCounter9TimesUs(i))+500000);
        sprintf('filled an empty itiCounterTime')
    end
end
%%
cerebellarFrameAnalyzer;

%dataStruct.cameraFrameTimes = get_frame_time_by_movie_info(ds)
%%
sum(dataStruct.puffRun)
sum(dataStruct.puffStill)
dataStruct.image = readtiff(pathName);
%% mean image field things
dataStruct.avgF = squeeze(mean(mean(dataStruct.image,1),2))';
%% end things
dataStruct.behaviorFramesAsCameraFrames = cumsum(dataStruct.goodFramesIx); %each space represents which camera frame it is. Each value represents which 
frameWindow = 15;
stimOns = find(dataStruct.stimOnIx);
%plot average fluoresence of each timepoint in a trial
for i = 1:length(stimOns);
    stim = stimOns(i);
    frameBeg = stim-frameWindow;   
    frameEnd = stim+frameWindow;
    frames(i,:) = dataStruct.avgF(frameBeg:frameEnd);
end
dataStruct.stimOnClips = frames;
dataStruct.stimPlot=mean(dataStruct.stimOnClips,1);
%creates a df/f movie to later average
Fn = squeeze(mean(dataStruct.image(:,:,[400:500]),3));     %HARDCODED  enter in an identified quiescent period into the 3rd dim of dataStruct.image
dFoverF = NaN(size(dataStruct.image));
%for i = 45:size(dFoverF,3)
for i = 1:size(dFoverF,3);
    dFoverF(:,:,i) = (double(dataStruct.image(:,:,i))-Fn)./Fn;
end
dataStruct.dFoverF = dFoverF;
%dataStruct.dFoverF = dFoverF(45:end);

%create a movie showing avgF before during and after puff
avgMovie = NaN(1002,1004,[2*frameWindow+1],length(stimOns));
aa = 1; 
for i = stimOns;
    avgMovie(:,:,:,aa) = dataStruct.image(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
avgMovie = squeeze(mean(avgMovie,4));
%wrtie a tiff file to view avgMovie    
%writetiff(avgMovie,[pathName '\avgMovie']);    
    
%create avg trace of the dfoverf
dFoverFAvg = NaN(1002,1004,[2*frameWindow+1],length(stimOns));
aa = 1; 
for i = stimOns;
%for i = stimOns(2:end)-44;
    dFoverFAvg(:,:,:,aa) = dataStruct.dFoverF(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
dFoverFAvg = squeeze(mean(dFoverFAvg,4));
%wrtie a tiff file to view    
%writetiff(dFoverFAvg,[pathName '\dFoverFAvg']);  

%%puffRun and puffStill averages
runPuffFrameNums = find(dataStruct.puffRun);
avgPuffRun = NaN(1002,1004,[2*frameWindow+1],length(runPuffFrameNums));
aa = 1; 
for i = runPuffFrameNums;
    avgPuffRun(:,:,:,aa) = dataStruct.image(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
avgPuffRun = squeeze(mean(avgPuffRun,4));
writetiff(avgPuffRun, [pathName '\avgPuffRun']);

stillPuffFrameNums = find(dataStruct.puffStill);
avgPuffStill = NaN(1002,1004,[2*frameWindow+1],length(stillPuffFrameNums));
aa = 1; 
for i = stillPuffFrameNums
    avgPuffStill(:,:,:,aa) = dataStruct.image(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
avgPuffStill = squeeze(mean(avgPuffStill,4));
writetiff(avgPuffStill, [pathName '\avgPuffStill']);    


%still minus run
stillMinusRun = avgPuffStill - avgPuffRun;
writetiff(stillMinusRun, [pathName '\stillMinusRun']);

%%plotting running vs F
B = [1:length(dataStruct.goodFramesIx)-2];
figure; plotyy(B,4*dataStruct.locomotionMatFrames(3:end),B,dataStruct.avgF(3:length(B)+2))
sum(puffRun)
sum(puffStill)