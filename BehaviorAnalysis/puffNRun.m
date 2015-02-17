% function [dataStruct] = puffNRun(folder)
imgName = 'Z:\Data\WidefieldImaging\GCaMP\150202_img19_2\150202_img19_2_MMStack.ome';
input = load('Z:\Data\WidefieldImaging\GCaMP\150202_img19_2\data-img19-150203-1714');
if isfield(input,'input'),
    input = input.input;
end
info = imfinfo(imgName, 'TIF');
cerebellarFrameAnalyzer
%dataStruct.cameraFrameTimes = get_frame_time_by_movie_info(ds)
%%
pathName = 'Z:\Data\WidefieldImaging\GCaMP\150202_img19_2';
dataStruct.image = readtiff(pathName);
%% mean image field things
dataStruct.avgF = squeeze(mean(mean(dataStruct.image,1),2))';
%% end things
dataStruct.behaviorFramesAsCameraFrames = cumsum(dataStruct.goodFramesIx);
frameWindow = 10;
stimOns = find(dataStruct.stimOnIx);
%plot average fluoresence of each timepoint in a trial
for i = 1:length(stimOns),
    stim = stimOns(i);
    frameBeg = stim-frameWindow;
    frameEnd = stim+frameWindow;
    frames(i,:) = dataStruct.avgF(frameBeg:frameEnd);
end
dataStruct.stimOnClips = frames;
dataStruct.stimPlot=mean(dataStruct.stimOnClips,1);
%creates a df/f movie to later average
Fn = squeeze(mean(dataStruct.image,3));
dFoverF = NaN(size(dataStruct.image));
for i = 1:size(dFoverF,3);
    dFoverF(:,:,i) = (double(dataStruct.image(:,:,i))-Fn)./Fn;
end
dataStruct.dFoverF = dFoverF;

%create a movie showing avgF before during and after puff
avgMovie = NaN(1002,1004,[2*frameWindow+1],length(stimOns));
aa = 1; 
for i = stimOns;
    avgMovie(:,:,:,aa) = dataStruct.image(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
avgMovie = squeeze(mean(avgMovie,4));
%wrtie a tiff file to view avgMovie    
writetiff(avgMovie,'Z:\Data\WidefieldImaging\GCaMP\150202_img19_2\avgMovie');    
    
%create avg trace of the dfoverf
dFoverFAvg = NaN(1002,1004,[2*frameWindow+1],length(stimOns));
aa = 1; 
for i = stimOns;
    dFoverFAvg(:,:,:,aa) = dataStruct.dFoverF(:,:,[i-frameWindow:i+frameWindow]);
    aa = aa +1;
end
dFoverFAvg = squeeze(mean(dFoverFAvg,4));
%wrtie a tiff file to view    
writetiff(dFoverFAvg,'Z:\Data\WidefieldImaging\GCaMP\150202_img19_2\dFoverFAvg');  
    
