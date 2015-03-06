% function [dataStruct] = puffNRun(folder)
pathName = 'Z:\Data\WidefieldImaging\GCaMP\150218_img15_1';
imgName = [pathName '\150218_img15_1_MMStack.ome'];
input = load([pathName '\data-img15-150218-1556']);
if isfield(input,'input'),
    input = input.input;
end
info = imfinfo(imgName, 'TIF');
cerebellarFrameAnalyzer;
%dataStruct.cameraFrameTimes = get_frame_time_by_movie_info(ds)
%%
%sum(dataStruct.puffRun)
%sum(dataStruct.puffStill)
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
Fn = squeeze(mean(dataStruct.image(:,:,[600:700]),3));     %HARDCODED  enter in an identified quiescent period into the 3rd dim of dataStruct.image
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

