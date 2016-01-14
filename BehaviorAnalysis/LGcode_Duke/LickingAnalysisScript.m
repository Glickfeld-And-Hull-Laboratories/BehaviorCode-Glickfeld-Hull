clear all
data_base = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
out_base = 'Z:\home\lindsey\Analysis\Behavior\Court\Licking_data';
date_mat = strvcat('150716','150717','150719');
mouse = '927';
for idate = 1:size(date_mat,1)
    date = date_mat(idate,:);
    dest = fullfile(out_base,[date '_' mouse]);
    data_dest = fullfile(data_base,['data-i' mouse '-' date '-*']);
    n = dir(data_dest);
    load(fullfile(data_base,n.name))

    ds = input;
    ntrials =size(cell2mat(ds.holdTimesMs),2);
    trialDurationMs = zeros(1,ntrials);
    holdStartsMs = zeros(1,ntrials);
    holdEndsMs = zeros(1,ntrials);
    postHoldMs = zeros(1,ntrials);
    lickTimesUs = ds.lickometerTimesUs;
    lickVals = ds.lickometerValues;

    for trial = 1:ntrials
        if size(lickTimesUs{trial},2)>0
            trialDurationMs(:,trial) =  (lickTimesUs{trial}(end)/1000) - (lickTimesUs{trial}(1)/1000);
            holdStartsMs(:,trial) = ds.holdStartsMs{trial} - (lickTimesUs{trial}(1)/1000);
            holdEndsMs(:,trial) = ds.holdStartsMs{trial} - (lickTimesUs{trial}(1)/1000) + ds.holdTimesMs{trial};
            postHoldMs(:,trial) = (lickTimesUs{trial}(end)/1000) - ds.holdStartsMs{trial};
        else
            trialDurationMs(:,trial) =  NaN;
            holdStartsMs(:,trial) = NaN;
            holdEndsMs(:,trial) = NaN;
            postHoldMs(:,trial) = NaN;
        end
    end
    max_trial = max(holdStartsMs,[],2)+max(postHoldMs,[],2)+1;
    max_start = max(holdStartsMs,[],2)+1;
    max_end = max(holdEndsMs,[],2)+1;
    lever_tc_press = zeros(max_trial, ntrials);
    lever_tc_release = zeros(max_trial, ntrials);

    for trial = 1:ntrials
        if size(lickTimesUs{trial},2)>0 
            lickTimeMs_press = (lickTimesUs{trial} - lickTimesUs{trial}(1))/1000 - holdStartsMs(:,trial)+max_start;
            lickTimeMs_release = (lickTimesUs{trial} - lickTimesUs{trial}(1))/1000 - holdEndsMs(:,trial)+max_end;
            lever_tc_press(lickTimeMs_press,trial) = 1;
            lever_tc_release(lickTimeMs_release,trial) = 1; 
        end
    end

    missedIx = strcmp(ds.trialOutcomeCell, 'ignore');
    failIx = strcmp(ds.trialOutcomeCell, 'failure');
    successIx = strcmp(ds.trialOutcomeCell, 'success');

%     figure
%     tr = [1:max_trial]-max_end;
%     for trial = 1:ntrials
%         plot(tr(1,max_end-1000:max_end+1000), lever_tc_release(max_end-1000:max_end+1000,trial), '-k')
%         hold on
%     end
%     title('All Releases')
%     xlabel('Time (ms)')
%     ylabel('Licks')
%     print([dest '_all_release_overlay.eps'], '-depsc');
%     print([dest '_all_release_overlay.pdf'], '-dpdf');
%     savefig([dest '_all_release_overlay.fig']);
    tr = [1:max_trial]-max_end;
    win = tr(1,max_end-1000:max_end+999);
    win_down = squeeze(mean(reshape(win, [50 40]),1));
    lever_release_win = lever_tc_release(max_end-1000:max_end+999,:);
    lever_release_down = squeeze(mean(reshape(lever_release_win, [50 40 ntrials]),1))*1000;
    figure
    avg_release = nanmean(lever_release_down,2);
    sem_release = nanstd(lever_release_down,[],2)./sqrt(sum(~isnan(lever_release_down(1,:)),2));
    shadedErrorBar(win_down, avg_release, sem_release, '-k');
    title([date ' ' mouse ' - All Releases- avg +/- sem - n =' num2str(ntrials)])
    xlabel('Time (ms)')
    ylabel('Lick rate (Hz)')
    print([dest '_release_avg.eps'], '-depsc');
    print([dest '_release_avg.pdf'], '-dpdf');
    savefig([dest '_release_avg.fig']);

    figure
    avg_success = nanmean(lever_release_down(:,find(successIx)),2);
    sem_success = nanstd(lever_release_down(:,find(successIx)),[],2)./sqrt(sum(~isnan(lever_release_down(:,find(successIx))),2));
    shadedErrorBar(win_down, avg_success, sem_success, '-k');
    hold on
    avg_fail = nanmean(lever_release_down(:,find(failIx)),2);
    sem_fail = nanstd(lever_release_down(:,find(failIx)),[],2)./sqrt(sum(~isnan(lever_release_down(:,find(failIx))),2));
    shadedErrorBar(win_down, avg_fail, sem_fail, '-r');
    title([date ' ' mouse ' - Successes- Black- n =' num2str(sum(successIx)) '; Failures- Red- n = ' num2str(sum(failIx))])
    xlabel('Time (ms)')
    ylabel('Lick rate (Hz)')
    print([dest '_success_fail_avg.eps'], '-depsc');
    print([dest '_success_fail_avg.pdf'], '-dpdf');
    savefig([dest '_success_fail_avg.fig']);

    save([dest '_lick.mat'], 'lever_tc_release', 'lever_release_down', 'win', 'win_down', 'max_start', 'max_end', 'avg_release','avg_success', 'avg_fail');
end