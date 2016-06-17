%Script for reading data from Running v. Fluorescence graphs produced by
%newpuffNRun. This script also calculates the average dF/F values for
%running and still states.

h=gcf;
axesObjs=get(h,'Children');
dataObjs=get(axesObjs,'Children');
dFoverbaseF = get(dataObjs{1},'YData');
%dF = dF{1};
runSpeed = get(dataObjs{2},'YData');
runSpeed = runSpeed{end};
rundF = dF(abs(runSpeed)>=4);
stilldF = dF(abs(runSpeed)<4);
rundF= mean(rundF);
stilldF = mean(stilldF);

disp(['Average dF/F while running: ' num2str(rundF)]);
disp(['Average dF/F while still: ' num2str(stilldF)]);


%% Find correlation between running speeds and dF/F
clear corrRunDf
h=gcf;
axesObjs=get(h,'Children');
dataObjs=get(axesObjs,'Children');
dF = get(dataObjs{2},'YData');
dF = cell2mat(dF); 
dF = flipud(dF);
runSpeed = get(dataObjs{3},'YData');
runSpeed = runSpeed{end};

figure;

dF = totaldF;
unqRunSpeeds = unique(totalRS); %unqRunSpeeds(end-1:end)=[];
colors = {'g', 'b', 'm', 'k', 'r', 'y', 'c', 'g','b'};
corrRunDf = zeros(length(unqRunSpeeds),3);
for jj=1:size(dF,1)
    for ii=1:length(unqRunSpeeds)
    runSpeedInd = find(totalRS == unqRunSpeeds(ii));
    corrRunDf(ii,1) = unqRunSpeeds(ii);
    corrRunDf(ii,2) = mean(totaldF(jj,runSpeedInd));
    corrRunDf(ii,3) = std(totaldF(jj,runSpeedInd))/length(runSpeedInd);
    end
    errorbar(corrRunDf(:,1),corrRunDf(:,2),corrRunDf(:,3),colors{jj})
    hold on
end

title('dF/F vs. Running Speed')
xlabel('Running Speed (cm/s)')
ylabel('dF/F')
%legend({'ROI 1', 'ROI 2', 'ROI 3', 'ROI 4', 'ROI 5'})
xlim([min(corrRunDf(:,1))-1 max(corrRunDf(:,1))+1])

%Firing rate vs run speed tuning curve
figure;
corrRunDf = zeros(length(unqRunSpeeds),3);
for jj=1:size(allfrate,1)
    for ii=1:length(unqRunSpeeds)
    runSpeedInd = find(runSpeed == unqRunSpeeds(ii));
    corrRunDf(ii,1) = unqRunSpeeds(ii);
    corrRunDf(ii,2) = mean(allfrate(jj,runSpeedInd));
    corrRunDf(ii,3) = std(allfrate(jj,runSpeedInd))/length(runSpeedInd);
    end
    errorbar(corrRunDf(:,1),corrRunDf(:,2),corrRunDf(:,3))
    hold on
end


runStim = zeros(1,length(runMat));
stimRunMat = runMat; stimRunMat(logical(dataStruct.stimStartStopIx))=2; %creates a vector of frames labelled with a 1 for run or still, and a 2 for the start/stop of stimulus

totalStim = [];
totalRS = [];
totaldF = [];


stimRunMat = zeros(1,length(stimOnIx));
stimRunMat(strfind(stimOnIx,[0 1])+1) = 2;
stimRunMat(strfind(stimOnIx,[1 0])) = 2;
totalStim = [totalStim stimRunMat];
totalRS = [totalRS runSpeed];
totaldF = [totaldF dFoverbaseF];


%stimulus-triggered tuning curve
stimFrames = find(stimRunMat==2);
brakedF = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
brakeRun = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
xx = 1;
goodstimind = find(diff(stimFrames)==stimFrames(2)-stimFrames(1));
for ii=1:length(find(diff(stimFrames)==stimFrames(2)-stimFrames(1)))
    brakedF(xx,:) = dFoverbaseF(1,stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    brakeRun(xx,:) = runSpeed(stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    xx=xx+1;
end

figure;
unqBrakeRun = unique(brakeRun);
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end

errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'r')

brakedF = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
brakeRun = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
xx = 1;
for ii=1:length(find(diff(stimFrames)==stimFrames(2)-stimFrames(1)))
    brakedF(xx,:) = dFoverbaseF(2,stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    brakeRun(xx,:) = runSpeed(stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    xx=xx+1;
end
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end
errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'k')

%% 
randomframes = randi(length(totaldF),1,400);
randomframes2 = randomframes +10;
randomframes3 = zeros(1,2*length(randomframes));
yy=1;
for ii=1:2:length(randomframes3)
    randomframes3(ii) = randomframes(yy);
    randomframes3(ii+1) = randomframes2(yy);
    yy=yy+1;
end

stimFrames = find(totalStim==2);
stimFrames = randomframes3;

brakedF = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1-25);
brakeRun = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1-25);
xx = 1;
for ii=1:2:length(stimFrames)
    brakedF(xx,:) = totaldF(1,stimFrames(ii):stimFrames(ii+1));
    brakeRun(xx,:) = totalRS(stimFrames(ii):stimFrames(ii+1));
    xx=xx+1;
end


unqBrakeRun = unique(brakeRun);
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end

errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'m')

brakedF = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1-25);
brakeRun = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1-25);
xx = 1;
for ii=1:2:length(stimFrames)
    brakedF(xx,:) = totaldF(2,stimFrames(ii):stimFrames(ii+1));
    brakeRun(xx,:) = totalRS(stimFrames(ii):stimFrames(ii+1));
    xx=xx+1;
end
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end
errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'c')
%% 
xx = 1;
goodstimind = find(diff(stimFrames)==stimFrames(2)-stimFrames(1));
for ii=1:length(find(diff(stimFrames)==stimFrames(2)-stimFrames(1)))
    brakedF(xx,:) = totaldF(1,stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    brakeRun(xx,:) = totalRS(stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    xx=xx+1;
end

unqBrakeRun = unique(brakeRun);
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end

errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'r')

brakedF = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
brakeRun = zeros(length(stimFrames)/2,stimFrames(2)-stimFrames(1)+1);
xx = 1;
for ii=1:length(find(diff(stimFrames)==stimFrames(2)-stimFrames(1)))
    brakedF(xx,:) = totaldF(2,stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    brakeRun(xx,:) = totalRS(stimFrames(goodstimind(ii)):stimFrames(goodstimind(ii)+1));
    xx=xx+1;
end
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end
errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'k')

%%
stimFrames = find(totalStim==2);
stimFrames(end)=[];
brakedF = zeros(length(stimFrames),15);
brakeRun = zeros(length(stimFrames),15);
xx = 1;
for ii=1:length(stimFrames)
    brakedF(xx,:) = totaldF(1,stimFrames(ii):stimFrames(ii)+14);
    brakeRun(xx,:) = totalRS(stimFrames(ii):stimFrames(ii)+14);
    xx=xx+1;
end


unqBrakeRun = unique(brakeRun);
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end

errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'r')

brakedF = zeros(length(stimFrames),15);
brakeRun = zeros(length(stimFrames),15);
xx = 1;
for ii=1:length(stimFrames)
    brakedF(xx,:) = totaldF(2,stimFrames(ii):stimFrames(ii)+14);
    brakeRun(xx,:) = totalRS(stimFrames(ii):stimFrames(ii)+14);
    xx=xx+1;
end
for jj=1:length(unqBrakeRun)
    brakeRunInd = find(brakeRun==unqBrakeRun(jj));
    brakeCorr(jj,1) = unqBrakeRun(jj);
    brakeCorr(jj,2) = mean(brakedF(brakeRunInd));
    brakeCorr(jj,3) = std(brakedF(brakeRunInd))/length(brakeRunInd);
end
errorbar(brakeCorr(:,1),brakeCorr(:,2),brakeCorr(:,3),'k')
