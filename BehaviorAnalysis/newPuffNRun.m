%Old puffNRun adapted for new cerebellarStim program.
% function [dataStruct] = puffNRun(folder)
clear; 
pathName = 'S:\VermisImaging\160219_img40_3\';
imgName = [pathName '160219_img40_3_MMStack.ome'];
input1 = load([pathName 'data-i940-160219-1705']);
if isfield(input1,'input'),
    input1 = input1.input;
end
info = imfinfo(imgName, 'TIF');

%%
tic;
newCerebellarFrameAnalyzer;
CFAtime = toc;
tic;
roi_selector;
ROItime = toc;
%%
dataStruct.image = readtiff([imgName '.tif']);
tic;
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
maskmovietime = toc;    
%% mean image field things
dataStruct.avgF = squeeze(mean(mean(dataStruct.image,1),2))';
%% end things
tic;
dataStruct.behaviorFramesAsCameraFrames = cumsum(dataStruct.goodFramesIx); %each space represents which camera frame it is. Each value represents which 
frameWindow = 1;
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

tifftime = toc;

%% Calculate base F
tic;
%Uses runMat matrix to return an index which adds up consecutive 1s. ie
%[0 1 1 0 1 1 1] = [0 1 2 0 1 2 3]
runMat(8000:8040) = 1;
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
if max(consecRunIx) > 30
    for ii=1:cluster.num_cluster
    maskBaseF{ii}= mean(dataStruct.(sprintf('avgF%d',ii))((longestRunIx(1)-26):longestRunIx(1)));
    end
end
if max(consecRunIx) == 30
    for ii=1:cluster.num_cluster
    maskBaseF{ii}= mean(dataStruct.(sprintf('avgF%d',ii))((longestRunIx(1)-28):longestRunIx(1)-2));
    end
end
%Ends program if there isn't at least a 3 second period of locomotion.
if max(consecRunIx) < 30
    disp('No 3 second period of consecutive locomotion found.');
    return
end

baseFcalctime = toc;
    
%% Plotting running vs F
tic;
B = [1:(length(dataStruct.goodFramesIx)-2)]; %cut out first two locomotion mat frames because of start artifact. Therefore need to cut out first two imaging frames even though they are good frames
B = B*(frameIntMs/1000);
runvFfig = figure;
%Set/calculate variables with proper units
revoLength = 30; revoPulse = 128;    %revoLength = length of 1 revolution (cm). revoPulse = # of pulses/revolution.
runSpeed = (1000/frameIntMs)*(revoLength/revoPulse)*dataStruct.locomotionMatFrames(3:end); %runSpeed=(frames/sec)*(cm/pulse)*(pulses/frame)
baseF = mean(dataStruct.avgF(3:length(B)+2));
for ii=1:cluster.num_cluster
    dFoverbaseF(ii,:) = (dataStruct.(sprintf('avgF%d',ii))(3:length(B)+2)/maskBaseF{ii})-1;
    B2(ii,:) = B;
end
%Plot Running v Fluorescence graph with appropriately labeled axes
[a b c] = plotyy(B,runSpeed,B2',dFoverbaseF')
title('Running vs. Fluorescence')
set(a(1), 'Ylim', [0 (max(runSpeed)+3)])   %runSpeed y-axis from 0 to max runSpeed (+3 to look nicer) 
set(a(2), 'Ylim', [min(min(dFoverbaseF)) max(max(dFoverbaseF))])    %dF/F y-axis from min dF/F to max dF/F
set(a(1), 'XTick', [0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]) %Converts x-axis frame ticks to seconds
set(a(2), 'XTick', [0,100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
set(a(1), 'XLim', [0 max(B)])  %Sets x-axis limits from 0 to end of data
set(a(2), 'XLim', [0 max(B)])
ylabel('Running Speed(cm/s)')
xlabel('Time (s)')
set(get(a(2),'Ylabel'),'String','dF/F')


%Plot dF/F over time for stimuli during Running
rundFvtimefig = figure;
stimWindow = 3000/frameIntMs; %Calculates # of frames to equal 3 seconds (span of stimulus window) based on frame interval.
stimTime = [0:2*stimWindow];    %x-axis for plot: stimulus window
meanrunStimdFoverF = NaN(length(runPuffFrameNums)-1,2*stimWindow+1,cluster.num_cluster);
for ii=1:(length(runPuffFrameNums)-1)
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
meanstillStimdFoverF = NaN(length(stillPuffFrameNums)-1,2*stimWindow+1,cluster.num_cluster);
for ii=(1:length(stillPuffFrameNums)-1)
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

graphingtime = toc;

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
save([pathName 'fulldata.mat'],'locoMat','blankIx','cluster','dFoverbaseF','runMat','runSpeed','stillMat','stimOnIx')

