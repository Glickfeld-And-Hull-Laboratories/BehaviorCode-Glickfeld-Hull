%% choose mouse
mouse = '706';
nSessBlk = 3;
firstTrialsCutoff = 100;
%%
% open days of experiments list
rc = behavConstsAV;
fn = fullfile(rc.ashleyAnalysis,mouse,'behavior');
data_tbl = readtable(fullfile(fn,'training_timesheet.csv'));
data_tbl.date = num2str(data_tbl.date);
data_tbl.time = num2str(data_tbl.time);
nexp = size(data_tbl,1);


%%
doPlot = 0;
% antiLicks_4D_fig = figure;
% antiLicks_1D_fig = figure;
nLicks = zeros(4,nexp);
nLicks_err = zeros(4,nexp);
expt_4D = zeros(1,nexp);
for iexp = 1:nexp
    % load expt data
    expDate = data_tbl.date(iexp,:);
    expTime = data_tbl.time(iexp,:);
    title_str = [mouse '-' expDate];
    mw = ['data-i' title_str '-' expTime];
    load(fullfile(rc.behavData,mw))
    
    %%
    nTrials = double(input.trialSinceReset);
    nTrPerBlk = floor(nTrials/nSessBlk);
    stimOn = double(input.dGratingDurationMs);
    rewDelay = double(input.rewardDelayMs);
    rewWinMs = rewDelay + 4000; % length of window to count licks should be length of reward delay + 4s
    nHistBins_250ms = floor(rewWinMs/500);
    trDirection = cell2mat(input.tGratingDirectionDeg);
    [trDir_ind dirs] = findgroups(trDirection);
    rewDirs = cell2mat(input.rewardStim);
    rewDirs_ind = ismember(dirs,rewDirs);
    trRew = cell2mat(input.tRewardTrial);
    % training with 1 or 4 (2/2 rew/nonrew) directions?
    if length(dirs) > 1
        expt_4D(iexp) = logical(1);
    end
        
    lickTimes_stimOnMs = cellfun(@(x,y) (double(x)-double(y))/1000, input.lickometerTimesUs, input.stimOnUs, 'unif',0);

    nLicks_rewWin_ind = cellfun(@(x) x >= 0 & x < rewWinMs, lickTimes_stimOnMs, 'unif', 0);
    nLicks_rewWin = cellfun(@(x,y) length(x(y)), lickTimes_stimOnMs, nLicks_rewWin_ind);

    nLicks_delayWin_ind = cellfun(@(x) x >= 0 & x < rewDelay, lickTimes_stimOnMs, 'unif', 0);
    nLicks_delayWin = cellfun(@(x,y) length(x(y)), lickTimes_stimOnMs, nLicks_delayWin_ind);

    lickTimes_rewWin = cellfun(@(x,y) x(y), lickTimes_stimOnMs, nLicks_rewWin_ind, 'unif',0);

    blk_colors = brewermap(nSessBlk+2,'Blues');

    
    for idir = 1:length(dirs)
        hold on
        ind = trDir_ind == idir;
        if length(ind) > firstTrialsCutoff
            ind = ind(1:firstTrialsCutoff);
        end
        templick = nLicks_delayWin(ind);
        nLicks(idir,iexp) = mean(templick);
        nLicks_err(idir,iexp) = ste(templick,2);
    end
        
 if doPlot
    % average number of licks in delay window (stim on to reward delivery)
    setFigParams4Print('landscape');figure;
    suptitle({title_str;'anticipatory licks per block of trials'})
    for iblk = 1:nSessBlk
        tr_ind = (nTrPerBlk *(iblk-1)) +1: nTrPerBlk*iblk;
        nLicks = zeros(1,length(dirs));
        nLicks_err = zeros(1,length(dirs));
        for idir = 1:length(dirs)
            hold on
            ind = trDir_ind(tr_ind) == idir;
            templick = nLicks_delayWin(tr_ind);
            templick = templick(ind);
            nLicks(idir) = mean(templick);
            nLicks_err(idir) = ste(templick,2);
        end
        subplot(1,3,iblk)
        h = bar(dirs,nLicks);
        hold on
        errorbar(dirs,nLicks,nLicks_err,'k.');
        figXAxis(h.Parent,'stim direction',[-20 155],dirs,dirs)
        figYAxis(h.Parent,'n licks before reward',[0 2])
        figAxForm(h.Parent)
        title(['trials ' num2str(tr_ind(1)) '-' num2str(tr_ind(end))])
    end   

    figure;
    suptitle(title_str)
    for idir = 1:length(dirs)
            subplot(2,2,idir)
        for iblk = 1:nSessBlk
            tr_ind = (nTrPerBlk *(iblk-1)) +1: nTrPerBlk*iblk;
            hold on
            ind = trDir_ind(tr_ind) == idir;
            templick = lickTimes_rewWin(tr_ind);
            templick = cell2mat(templick(ind));
            h = cdfplot(templick);
            h.Color = blk_colors(iblk+2,:);
            h.LineWidth = 1;
        end
    figXAxis(h.Parent,'n licks in reward window',[0 rewWinMs])
    figYAxis(h.Parent,'fraction of trials',[])
    figAxForm(h.Parent)
    vline(rewDelay,'k--')
    title([num2str(dirs(idir)) ' deg'])
    leg = legend(strread(num2str(1:nSessBlk),'%s'),'location','southeast');
    title(leg,'session block')
    end
 end
end

figure;
subplot 121
h = bar(dirs,nLicks(:,logical(expt_4D)),'grouped');
hold on
% errorbar(dirs,nLicks(:,expt_4D),nLicks_err(:,expt_4D),'k.');
figXAxis([],'stim direction',[-20 155],dirs,dirs)
figYAxis([],'n licks before reward',[0 5])
figAxForm([])
title({'anticipitory licks'})

subplot 122
h = bar(nLicks(1,:));
hold on
errorbar(1:nexp,nLicks(1,:),nLicks_err(1,:),'k.');
figXAxis([],'training day',[0 nexp+1],1:nexp,1:nexp)
figYAxis([],'n licks before reward',[0 2])
figAxForm([])
title({'0 deg only'})

print(fullfile(fn,'antiLicks'),'-dpdf','-fillpage')