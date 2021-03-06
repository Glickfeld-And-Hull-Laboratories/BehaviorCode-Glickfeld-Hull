% function [dataStruct] = puffNRun(folder)
clear; 
pathName = 'S:\160113_img30_1\';
imgName = [pathName '160113_img30_1_MMStack.ome'];
input1 = load([pathName 'data-i930-160113-1417']);
if isfield(input1,'input'),
    input1 = input1.input;
end
info = imfinfo(imgName, 'TIF');

%This code is duplicated within cerebellarFrameAnalyzer
%%Fill any empty sp/itiCounters with a NaN
%Note: Sometimes stopAfterNTrials > length of sp/iti counters, suggesting
%that the experiment might have been stopped early. Good fix might be to
%set input.StopAfterNTrials = length of a counter

for i = 1:length(input1.itiCounter5)
     if isempty(cell2mat(input1.spCounter10(i)))==1
        input1.spCounter10(i) = input1.spCounter9(i);
        sprintf('filled an empty spCounter10')
    end
    if isempty(cell2mat(input1.itiCounter10(i)))==1
        input1.itiCounter10(i) = input1.itiCounter9(i);
        sprintf('filled an empty itiCounter10')
    end
end
%{
%Fill in any missing iti/spCounterTimes
for i = 1:length(input1.itiCounter5)
    if isempty(cell2mat(input1.spCounter10TimesUs(i)))
        input1.spCounter10TimesUs(i) =  mat2cell(cell2mat(input1.spCounter9TimesUs(i))+500000);
        sprintf('filled an empty spCounterTime')
    end
    if isempty(cell2mat(input1.itiCounter10TimesUs(i)))
        input1.itiCounter10TimesUs(i) =  mat2cell(cell2mat(input1.itiCounter9TimesUs(i))+500000);
        sprintf('filled an empty itiCounterTime')
    end
end
%}
%%
cerebellarFrameAnalyzer;
roi_selector;
%dataStruct.cameraFrameTimes = get_frame_time_by_movie_info(ds)
%%
dataStruct.image = readtiff([imgName '.tif']);
%For loop creates a "mask movie" for each ROI
for ii=1:cluster.num_cluster
    mask_index{ii} = NaN(1002,1004,size(dataStruct.image,3));
    for jj=1:size(dataStruct.image,3)
        mask_index{ii}(:,:,jj) = cluster.(sprintf('mask%d',ii));
    end
    mask_index{ii} = logical(mask_index{ii});
    mask_movie{ii} = dataStruct.image;
    mask_movie{ii}(~mask_index{ii})=0;
    dataStruct.(sprintf('avgF%d',ii))= squeeze(mean(mean(mask_movie{ii},1),2))';
end
    
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
Fn = squeeze(mean(dataStruct.image(:,:,[200:300]),3));     %HARDCODED  enter in an identified quiescent period into the 3rd dim of dataStruct.image
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

%% Calculate base F
%Uses runMat matrix to return an index which adds up consecutive 1s. ie
%[0 1 1 0 1 1 1] = [0 1 2 0 1 2 3]
mm = size(runMat,1);
consecRunIx = [reshape([zeros(1,mm);runMat.'],[],1);0];
pp=find(~consecRunIx); dd=1-diff(pp);
consecRunIx(pp)=[0;dd];
consecRunIx=reshape(cumsum(consecRunIx(1:end-1)),[],mm).';
consecRunIx(:,1)=[];
%Finds the location of the longest consecutive run.
longestRunIx= find(consecRunIx == max(consecRunIx));
%These if-statements make sure we get a buffered 3 second period of
%locomotion to calculate the baseF for each ROI, by averaging the avgF over
%those frames. 
if max(consecRunIx) > 12
    for ii=1:cluster.num_cluster
    maskBaseF{ii}= mean(dataStruct.(sprintf('avgF%d',ii))((longestRunIx(1)-26):longestRunIx(1)));
    end
end
if max(consecRunIx) == 12
    for ii=1:cluster.num_cluster
    maskBaseF{ii}= mean(dataStruct.(sprintf('avgF%d',ii))((longestRunIx(1)-28):longestRunIx(1)-2));
    end
end
%Ends program if there isn't at least a 3 second period of locomotion.
if max(consecRunIx) < 12
    disp('No 3 second period of consecutive locomotion found.');
    return
end

%% Plotting running vs F
B = [1:(length(dataStruct.goodFramesIx)-2)]; %cut out first two locomotion mat frames because of start artifact. Therefore need to cut out first two imaging frames even though they are good frames
B = B*(frameIntUs/1000000);
runvFfig = figure;
%Set/calculate variables with proper units
revoLength = 30; revoPulse = 32;    %revoLength = length of 1 revolution (cm). revoPulse = # of pulses/revolution.
runSpeed = (1000000/frameIntUs)*(revoLength/revoPulse)*dataStruct.locomotionMatFrames(3:end); %runSpeed=(frames/sec)*(cm/pulse)*(pulses/frame)
baseF = mean(dataStruct.avgF(3:length(B)+2));
dFoverbaseF = (dataStruct.avgF1(3:length(B)+2)/maskBaseF{1})-1;
%Plot Running v Fluorescence graph with appropriately labeled axes
[a b c] = plotyy(B,runSpeed,B,dFoverbaseF)
title('Running vs. Fluorescence')
set(a(1), 'Ylim', [0 (max(runSpeed)+3)])   %runSpeed y-axis from 0 to max runSpeed (+3 to look nicer) 
set(a(2), 'Ylim', [min(dFoverbaseF) max(dFoverbaseF)])    %dF/F y-axis from min dF/F to max dF/F
set(a(1), 'XTick', [0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]) %Converts x-axis frame ticks to seconds
set(a(2), 'XTick', [0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
set(a(1), 'XLim', [0 max(B)])  %Sets x-axis limits from 0 to end of data
set(a(2), 'XLim', [0 max(B)])
ylabel('Running Speed(cm/s)')
xlabel('Time (s)')
set(get(a(2),'Ylabel'),'String','dF/F')


%Plot dF/F over time for stimuli during Running
rundFvtimefig = figure;
stimWindow = 3000000/frameIntUs; %Calculates # of frames to equal 3 seconds (span of stimulus window) based on frame interval.
stimTime = [0:2*stimWindow];    %x-axis for plot: stimulus window
meanrunStimdFoverF = NaN(length(runPuffFrameNums),2*stimWindow+1,cluster.num_cluster);
for ii=1:length(runPuffFrameNums)
    runStimFrame = runPuffFrameNums(ii);
    for jj=1:cluster.num_cluster
    runStimdFoverF = ((dataStruct.(sprintf('avgF%d',jj))(runStimFrame-stimWindow:runStimFrame+stimWindow)/maskBaseF{jj})-1); %Gather F data for frames in desired window around stimulus
    meanrunStimdFoverF(ii,:,jj) = runStimdFoverF; %Build matrix of F data around stimulus for all trials
    end
end

%For loop averages dF over F for each ROI
for ii = 1:cluster.num_cluster
    roidFoverFrun = mean(meanrunStimdFoverF,1);
    plot(stimTime,roidFoverFrun(:,:,ii))
    hold on;
end

meanrunStimdFoverF = mean(mean(meanrunStimdFoverF,1),3); %Average dF/F over all trials and add to plot as red line
plot(stimTime,meanrunStimdFoverF,'r')
title('dF/F for Running Stimuli')
xlabel('Time (s)')
ylabel('dF/F')
set(gca,'XLim',[0 2*stimWindow]);
set(gca,'XTick', [0 stimWindow 2*stimWindow]);
set(gca,'XTickLabel',[-3 0 3]);

%Plot dF/F over time for stimuli while Still
stilldFvtimefig = figure;
meanstillStimdFoverF = NaN(length(stillPuffFrameNums),2*stimWindow+1,cluster.num_cluster);
for ii=1:length(stillPuffFrameNums)
    stillStimFrame = stillPuffFrameNums(ii);
    for jj=1:cluster.num_cluster
    stillStimdFoverF = ((dataStruct.(sprintf('avgF%d',jj))(stillStimFrame-stimWindow:stillStimFrame+stimWindow)/maskBaseF{jj})-1);
    meanstillStimdFoverF(ii,:,jj) = stillStimdFoverF;
    end
end

%For loop averages dF over F for each ROI
for ii=1:cluster.num_cluster
    roidFoverFstill = mean(meanstillStimdFoverF,1);
    plot(stimTime,roidFoverFstill(:,:,ii))
    hold on;
end

meanstillStimdFoverF = mean(mean(meanstillStimdFoverF,1),3);
plot(stimTime,meanstillStimdFoverF,'r')
title('dF/F for Still Stimuli')
xlabel('Time (s)')
ylabel('dF/F')
set(gca,'XLim',[0 2*stimWindow]);
set(gca,'XTick', [0 stimWindow 2*stimWindow]);
set(gca,'XTickLabel',[-3 0 3]);

%% Outputs number of trials w/ puffs delivered while running/still
disp(['Number of puffs while running: ' num2str(sum(puffRun))])
disp(['Number of puffs while still: ' num2str(sum(puffStill))])
disp([pathName])
disp([imagingRate])
disp([length(dataStruct.goodFramesIx)])

%% Save important figures automatically

saveas(runvFfig,[pathName 'RunvF'], 'fig')
saveas(rundFvtimefig,[pathName 'RunStim'], 'fig')
saveas(stilldFvtimefig,[pathName 'StillStim'], 'fig')








