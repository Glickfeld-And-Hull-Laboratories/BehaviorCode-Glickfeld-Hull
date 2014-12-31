dateIx314 = [27 29 31 32 36];
files314 = dataStru(314).check.downloadedname';

timingIx314 = dataStru(314).values.timing;
timingIx314 = timingIx314(dateIx314);
fCue314 = dataStru(314).values.fixedCue;
fCue314 = fCue314(dateIx314)
rCue314 = dataStru(314).values.randCue;
rCue314 = rCue314(dateIx314)
cueIx314 = fCue314|rCue314;

tStarts =lidoStruct(314).values.lidoStarts;
tROR = lidoStruct(314).values.reactOverRequired;
cueROR1 = tROR(cueIx314);
timingROR1 = tROR(logical(timingIx314));
cueStarts1 = tStarts(cueIx314);
timingStarts1 = tStarts(logical(timingIx314));


dateIx315 = [42 45 47 50];
files315 = dataStru(315).check.downloadedname'
 
timingIx315 = dataStru(315).values.timing;
timingIx315 = timingIx315(dateIx315);
fCue315 = dataStru(315).values.fixedCue;
fCue315 = fCue315(dateIx315)
rCue315 = dataStru(315).values.randCue;
rCue315 = rCue315(dateIx315)
cueIx315 = fCue315|rCue315;

tStarts =lidoStruct(315).values.lidoStarts;
tROR = lidoStruct(315).values.reactOverRequired;
cueROR2 = tROR(cueIx315);
timingROR2 = tROR(logical(timingIx315));
cueStarts2 = tStarts(cueIx315);
timingStarts2 = tStarts(logical(timingIx315));

cueStarts = [cell2mat(cueStarts1) cell2mat(cueStarts2)];
cueROR = [cell2mat(cueROR1) cell2mat(cueROR2)];
timingStarts = [cell2mat(timingStarts1) cell2mat(timingStarts2)];
timingROR = [cell2mat(timingROR1) cell2mat(timingROR2)];

[proCueStarts, I] = sort(cueStarts);
proCueROR = cueROR(I);



[proTimingStarts, I] = sort(timingStarts);
proTimingROR = timingROR(I);

[procStarts, I] = sort(rawStarts);
rawROR = cell2mat(lidoStruct(j).values.reactOverRequired);
procROR = rawROR(I);

plot(proCueStarts, smooth(proCueROR, 30, 'moving'))
 xlabel('Time from Lidocaine Onset (minutes)')
    ylabel('Reaction Time / Required Hold Time for Reward')
    tName = strcat('Lidocaine Timing Analysis - Averaged Trial Data - Cued Conditions');
    vH = vert_lines(0);
        set(vH, 'Color', 'g');
        set(vH, 'LineWidth', 2);
    title(tName)
    sName = strcat('~/Documents/MWorks/DailyPlots/lidocaineTiming-cued-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 8], ...
             'PrintUI', false };
    exportfig_print(epParams{:});
    close

    
    
plot(proTimingStarts, smooth(proTimingROR, 30, 'moving'))
xlabel('Time from Lidocaine Onset (minutes)')
    ylabel('Reaction Time / Required Hold Time for Reward')
    tName = strcat('Lidocaine Timing Analysis - Averaged Trial Data - Timing Conditions');
    vH = vert_lines(0);
        set(vH, 'Color', 'g');
        set(vH, 'LineWidth', 2);
    title(tName)
    sName = strcat('~/Documents/MWorks/DailyPlots/lidocaineTiming-timing-', datestr(today, 'yymmdd'), '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 8], ...
             'PrintUI', false };
    exportfig_print(epParams{:});
    close
