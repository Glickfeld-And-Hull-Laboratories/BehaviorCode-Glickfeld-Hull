function [ dataStruct ] = plotHoldAndReactTimeCDFs( subjMat, endMat )
%When pointed to the dataPath folder of your choice and using the subjects
%contained in your subjMat matrix, should plot training CDF progression for
%single mice and average across multiple mice. It will produce hold time
%CDFs and reaction time CDFs, with training progression steps grouped in
%5-day intervals.
global dataStruct
addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
% Data Path Formatting
disp('Loading data from MWorks Data folder...')
dataPath = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\';
nSubjs= length(subjMat);
dataStruct=struct;

%% This section searches dataPath for files specific to each animal for the specified days, then loads specified parts of data into struData.
for j = subjMat;
    desNames = dir([dataPath 'data-i' num2str(j) '*']);
    nNames = length(desNames);
    for iN = 1:nNames;
        tName = fullfile(dataPath, desNames(iN).name);
        dataStruct(j).check.downloadedname{1,iN}= tName; 
        ds=mwLoadData(tName, 'max');
        
        nTrials = length(ds.trialOutcomeCell);
        dataStruct(j).values.nTrials{1,iN} = nTrials;
        %% trial outcome calculations
        corrIx = strcmp(ds.trialOutcomeCell, 'success');
        dataStruct(j).values.perCorrect{1,iN} = sum(corrIx)/nTrials;
        igIx = strcmp(ds.trialOutcomeCell, 'ignore');
        dataStruct(j).values.perIgnore{1,iN} = sum(igIx)/nTrials;
        incIx = strcmp(ds.trialOutcomeCell, 'incorrect');
        dataStruct(j).values.perIncorrect{1,iN} = sum(incIx)/nTrials;
        earlyIx = strcmp(ds.trialOutcomeCell, 'failure');
        dataStruct(j).values.perEarly{1,iN} = sum(earlyIx)/nTrials;
        a = mean(cell2mat(ds.tTotalReqHoldTimeMs));
        b = mean(cell2mat(ds.holdTimesMs));
        dataStruct(j).values.meanReqHold{1,iN} = a;
        dataStruct(j).values.meanActualHold{1,iN} = b;
        reacts = cell2mat(ds.reactTimesMs);
        dataStruct(j).values.meanReact{1,iN} = mean(reacts(corrIx));
        
        dataStruct(j).values.rawHolds{1, iN} = cell2mat(ds.holdTimesMs);
        dataStruct(j).values.rawReq{1, iN} = cell2mat(ds.tTotalReqHoldTimeMs);
        dataStruct(j).values.rawReact{1, iN} = cell2mat(ds.reactTimesMs);
        
        if a<500,
            dataStruct(j).values.trainingPhase{1,iN} = 1;
        end
        if and(a>=500, a<1500),
            dataStruct(j).values.trainingPhase{1,iN} = 2;
        end
        if a>=1500,
            dataStruct(j).values.trainingPhase{1,iN} = 3;
        end
        
    end
    disp(strcat('Subject Data from i', num2str(j), ' Loaded'))
end

beep
disp('Loading complete!')%%

%plot for number of trials by day
figure;
nTrials = cell(length(subjMat));
for j = 1:length(subjMat)
    sub = subjMat(j);
    nTrials{j} = cell2mat(dataStruct(sub).values.nTrials);
    nDays(j) = length(dataStruct(sub).values.nTrials);
    plot(nTrials{j},'ok')
    hold on
end
max_days = max(nDays,[],2);
trial_mat = NaN(length(subjMat),max_days);
for j = 1:length(subjMat)
    trial_mat(j,1:nDays(1,j)) = nTrials{j};
end
%plot minumum number of trials across mice for a day
for i = 1:max_days
    plot(i,min(trial_mat(:,i),[],1),'or')
    hold on
end


for j = 1:length(subjMat)
    sub = subjMat(j);
    desNames = dir([dataPath 'data-i' num2str(sub) '*']);
    nNames = size(desNames,1);
    dataStruct(sub).date = zeros(1,nNames);
    for iN = 1:nNames;
        dataStruct(sub).date(1,iN)= str2num(desNames(iN).name(11:16));
    end
    dataStruct(sub).endDay = find(dataStruct(sub).date == endMat(1,j));
end

trialMat = [100 400 400 400];
exMat = [0 50 50 50];
for j = 1:length(subjMat)
    sub = subjMat(j);
    dataStruct(sub).values.selectHolds= cell(1,4);
    dataStruct(sub).values.selectReact= cell(1,4);
    endDay(j) = dataStruct(sub).endDay;
    days = [1 10 19 endDay(j)];
    for id = 1:4;
        day = days(id);
        if dataStruct(sub).values.nTrials{1,day}-exMat(id)>trialMat(id)
            ind = randsample(dataStruct(sub).values.nTrials{1,day}-exMat(id),trialMat(id));
            rawHolds = dataStruct(sub).values.rawHolds{1,day};
            rawReact = dataStruct(sub).values.rawReact{1,day};
            dataStruct(sub).values.selectHolds{1,id} = rawHolds(1,ind);
            dataStruct(sub).values.selectReact{1,id} = rawReact(1,ind);
        else
            dataStruct(sub).values.selectHolds{1,id} = NaN;
            dataStruct(sub).values.selectReact{1,id} = NaN;
        end
        if id == 1
            maxhold(j) = max(rawHolds,[],2);
        end
    end
end

intervalHolds = cell(1,4);
intervalReacts = cell(1,4);
for ii = 1:4, % linspace(5:5:ceil(maxDays/5)*5),
    tHolds = zeros(1,trialMat(ii));
    tReact = zeros(1,trialMat(ii));
    for iii = 1:length(subjMat),
        sub = subjMat(iii);
        if ii==1
            if maxhold(iii) > 5000
                tHolds(iii,:) = dataStruct(sub).values.selectHolds{ii};
                tReact(iii,:) = dataStruct(sub).values.selectReact{ii};
            else
                tHolds(iii,:) = NaN(1,trialMat(ii));
                tReact(iii,:) = NaN(1,trialMat(ii));
            end
        else
            tHolds(iii,:) = dataStruct(sub).values.selectHolds{ii};
            tReact(iii,:) = dataStruct(sub).values.selectReact{ii};
        end
    end
    sza = size(tHolds);
    a = double(reshape(tHolds,[1, sza(1)*sza(2)]));
    intervalHolds{ii} = a;
    
    szb = size(tReact);
    b = double(reshape(tReact,[1, szb(1)*szb(2)]));
    intervalReacts{ii} = b;
end

colorCell = {[0 0 1], ... % blue
             [0 1 0], ... & green
             [0.8 0.8 0], ... % yellow
             [1 0.5 0], ... % orange
             [1 0 0], ... % red
             [0.8 0 0.8]}; % purple};
figure;
for i = 1:4,
    p1 = cdfplot(intervalReacts{i});
    set(p1, 'Color', colorCell{i});
    hold on
end
vert_lines(0, 'k')
xlim([-1000 5000])
ylim([0 1])
xlabel('Reaction Time (ms)')
title('Reaction Time CDFs')
legend('Day 1', 'Day 10', 'Day 20', ['Day ' num2str(chop(mean(endDay,2),2)) ' +/- ' num2str(chop(std(endDay,[],2)./sqrt(length(endDay)),2))],'Location', 'east')
savefig('Z:\home\lindsey\Analysis\Behavior\Court\HoldAndReactTimes\ReactSummaryCDF.fig')

figure;
for i = 1:4,
    p2 = cdfplot([intervalHolds{i}]);
    set(p2, 'Color', colorCell{i});
    hold on
end
xlim([0 5000])
ylim([0 1])
xlabel('Hold Time (ms)')
title('Hold Time CDFs')
legend('Day 1', 'Day 10', 'Day 20', ['Day ' num2str(chop(mean(endDay,2),2)) ' +/- ' num2str(chop(std(endDay,[],2)./sqrt(length(endDay)),2))],'Location', 'east')
savefig('Z:\home\lindsey\Analysis\Behavior\Court\HoldAndReactTimes\HoldSummaryCDF.fig')
end

