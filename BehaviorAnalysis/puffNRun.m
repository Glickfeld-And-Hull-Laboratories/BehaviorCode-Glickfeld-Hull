% function [dataStruct] = puffNRun(folder)
imgName = 'Z:\Data\WidefieldImaging\GCaMP\150212_img20_2\150212_img20_2_MMStack.ome';
input = load('Z:\Data\WidefieldImaging\GCaMP\150212_img20_2\data-img20-150212-1847.mat');
if isfield(input,'input'),
    input = input.input;
end
info = imfinfo(imgName, 'TIF');
%cerebellarFrameAnalyzer
%% EVENTUALLY this will compare the frame times from the camera and the behavior and make a harmonious dance
%dataStruct.cameraFrameTimes = get_frame_time_by_movie_info(ds)
%%
pathName = 'Z:\Data\WidefieldImaging\GCaMP\150212_img20_2';
dataStruct.image = readtiff(pathName);
dataStruct.behaviorFramesAsCameraFrames = cumsum(dataStruct.goodFramesIx);