addpath('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorAnalysis');
dataPath = '~/Documents/MWorks/Data';
dirListing = dir(dataPath);
fileNames= {dirListing.name};
fileNames = fileNames(5:end);
outR = regexp(fileNames, 'data-i([0-9]*)-([0-9]*)?', 'tokens');
outRV = cat(1, outR{:});
outRM = cat(1, outRV{:});
subjNumList = str2double(outRM(:,1));
dateStrList = str2double(outRM(:,2));

lidoStruct = struct;

lidocaineBook = '~/Desktop/LidocaineBook.xlsx';
[a, sheetNames] = xlsfinfo(lidocaineBook);
for i=1:length(sheetNames),
    j = str2num(sheetNames{i})
    subjIx = subjNumList == j;
    
    [num, txt, raw] = xlsread(lidocaineBook, sheetNames{i});
    dates = num(:,1);
    lidocaineOnsetTr = num(:,8);
    
    lidoStruct(j).values.lidocaineOnset = lidocaineOnsetTr;
    lidoStruct(j).values.dates = dates;
    
    for k=1:length(dates);
        m = dates(k);
        dateIx = dateStrList == m;
        subjDateIx = subjIx & dateIx;
        desNames = fileNames(subjDateIx);
        tName = fullfile(dataPath, desNames{1});
        ds=mwLoadData(tName, 'max');
        disp(ds.savedDataName)
        starts = double(cell2mat(ds.tThisTrialStartTimeMs));
        trialStartsMinutes = (starts - starts(1))/60000;
        lidoStruct(j).values.trialStartsMinutes{k} = trialStartsMinutes;
        lidoStarts = trialStartsMinutes - trialStartsMinutes(lidocaineOnsetTr(k));
        lidoStruct(j).values.lidoStarts{k} = lidoStarts;
        react = cell2mat(ds.reactTimesMs);
        req = cell2mat(ds.tTotalReqHoldTimeMs);
        reactOverRequired = double(react)./req;
        lidoStruct(j).values.reactOverRequired{k} = reactOverRequired;
        
        scatter(lidoStarts,reactOverRequired)
        vH = vert_lines(0);
        set(vH, 'Color', 'g');
        set(vH, 'LineWidth', 2);
        hold on
    end
    xlabel('Time from Lidocaine Onset (minutes)')
    ylabel('Reaction Time / Required Hold Time for Reward')
    tName = strcat('Lidocaine Timing Analysis - Raw Trial Data - i', num2str(j));
    title(tName)
    sName = strcat('~/Documents/MWorks/DailyPlots/lidocaineTiming-raw-i', num2str(j),'-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 8], ...
             'PrintUI', false };
    exportfig_print(epParams{:});
    
    close
    
    rawStarts = cell2mat(lidoStruct(j).values.lidoStarts);
    [procStarts, I] = sort(rawStarts);
    rawROR = cell2mat(lidoStruct(j).values.reactOverRequired);
    procROR = rawROR(I);
    
    plot(procStarts, smooth(procROR, 30, 'moving'))
    
    xlabel('Time from Lidocaine Onset (minutes)')
    ylabel('Reaction Time / Required Hold Time for Reward')
    tName = strcat('Lidocaine Timing Analysis - Averaged Trial Data - i', num2str(j));
    vH = vert_lines(0);
        set(vH, 'Color', 'g');
        set(vH, 'LineWidth', 2);
    title(tName)
    sName = strcat('~/Documents/MWorks/DailyPlots/lidocaineTiming-avg-i', num2str(j),'-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 8], ...
             'PrintUI', false };
    exportfig_print(epParams{:});
    close
end
%date = str2num(ds.saveTime(1:6));
%rw = find(dates==date);
%trialVals = raw(rw+1,:);