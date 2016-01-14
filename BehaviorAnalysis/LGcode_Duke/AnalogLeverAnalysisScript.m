clear all
data_base = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
out_base = 'Z:\home\lindsey\Analysis\Behavior\Court\Analog_lever';
date = '150708';
mouse = '318';
tick = 0.36;
dest = fullfile(out_base,[date '_' mouse '_']);
data_dest = fullfile(data_base,['data-i' mouse '-' date '-*']);
n = dir(data_dest);
load(fullfile(data_base,n.name))

ds = input;
ntrials =size(cell2mat(ds.holdTimesMs),2);
trialDurationMs = zeros(1,ntrials);
holdStartsMs = zeros(1,ntrials);
holdEndsMs = zeros(1,ntrials);
postHoldMs = zeros(1,ntrials);
quadTimesUs = ds.quadratureTimesUs;
quadVals = ds.quadratureValues;

for trial = 1:ntrials
    trialDurationMs(:,trial) =  (quadTimesUs{trial}(end)/1000) - (quadTimesUs{trial}(1)/1000);
    holdStartsMs(:,trial) = ds.holdStartsMs{trial} - (quadTimesUs{trial}(1)/1000);
    holdEndsMs(:,trial) = ds.holdStartsMs{trial} - (quadTimesUs{trial}(1)/1000) + ds.holdTimesMs{trial};
    postHoldMs(:,trial) = (quadTimesUs{trial}(end)/1000) - ds.holdStartsMs{trial};
end
max_trial = max(holdStartsMs,[],2)+max(postHoldMs,[],2)+1;
max_start = max(holdStartsMs,[],2)+1;
max_end = max(holdEndsMs,[],2)+1;
lever_tc_press = NaN(max_trial, ntrials);
lever_tc_release = NaN(max_trial, ntrials);

for trial = 1:ntrials
    quadTimeMs_press = (quadTimesUs{trial} - quadTimesUs{trial}(1))/1000 - holdStartsMs(:,trial)+max_start;
    quadTimeMs_release = (quadTimesUs{trial} - quadTimesUs{trial}(1))/1000 - holdEndsMs(:,trial)+max_end;
    quadVal = quadVals{trial};
    quadVal(quadVal>100) = NaN;
    lever_tc_press(quadTimeMs_press,trial) = quadVal;
    lever_tc_release(quadTimeMs_release,trial) = quadVal; 
end

lever_press_interp = NaN(max_trial, ntrials);
lever_release_interp = NaN(max_trial, ntrials);

for trial = 1:ntrials
    press_ind{trial} = find(~isnan(lever_tc_press(:,trial)));
    release_ind{trial} = find(~isnan(lever_tc_release(:,trial)));
    press_val{trial} = lever_tc_press(press_ind{trial},trial);
    release_val{trial} = lever_tc_release(release_ind{trial},trial);
    resamp_press = press_ind{trial}(1):press_ind{trial}(end);
    resamp_release = release_ind{trial}(1):release_ind{trial}(end);
    temp_press_interp = interp1(press_ind{trial},press_val{trial},resamp_press);
    temp_release_interp = interp1(release_ind{trial},release_val{trial},resamp_release);
    lever_press_interp(resamp_press,trial) = temp_press_interp;
    lever_release_interp(resamp_release,trial) = temp_release_interp;
end

lever_press_mm = lever_press_interp.*tick;
lever_release_mm = lever_release_interp.*tick;

good_start_trials = find(mean(lever_press_mm(max_start-150:max_start-50,:),1)>(-10*tick));
good_hold_trials = find(cell2mat(ds.holdTimesMs)>200);
good_press_trials = intersect(good_start_trials, good_hold_trials);

figure
tp = [1:max_trial]-max_start;
for trial = good_press_trials
    plot(tp(1,max_start-200:max_start+200), lever_press_mm(max_start-200:max_start+200,trial), '-k')
    hold on
end
title('Presses')
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_press_overlay.eps'], '-depsc');
print([dest '_press_overlay.pdf'], '-dpdf');
savefig([dest '_press_overlay.fig']);

avg_press = nanmean(lever_press_mm(max_start-200:max_start+200,good_press_trials),2);
sem_press = nanstd(lever_press_mm(max_start-200:max_start+200,good_press_trials),[],2)./sqrt(sum(~isnan(lever_press_mm(max_start-200:max_start+200,good_press_trials)),2));
figure; errorbar(tp(1,max_start-200:max_start+200)', avg_press, sem_press, '-k');
title(['Presses- avg +/- sem- n = ' num2str(length(good_press_trials))])
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_press_avg.eps'], '-depsc');
print([dest '_press_avg.pdf'], '-dpdf');
savefig([dest '_press_avg.fig']);

missedIx = strcmp(ds.trialOutcomeCell, 'ignore');
good_response_trials = find(missedIx==0);
good_end_trials = find(mean(lever_release_mm(max_end+100:max_end+200,:),1)>(-10*tick));
good_release_trials = intersect(intersect(good_hold_trials, good_response_trials),good_end_trials);

figure
tr = [1:max_trial]-max_end;
for trial = good_release_trials
    plot(tr(1,max_end-200:max_end+200), lever_release_mm(max_end-200:max_end+200,trial), '-k')
    hold on
end
title('Releases')
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_release_overlay.eps'], '-depsc');
print([dest '_release_overlay.pdf'], '-dpdf');
savefig([dest '_release_overlay.fig']);

avg_release = nanmean(lever_release_mm(max_end-200:max_end+200,good_release_trials),2);
sem_release = nanstd(lever_release_mm(max_end-200:max_end+200,good_release_trials),[],2)./sqrt(sum(~isnan(lever_release_mm(max_end-200:max_end+200,good_release_trials)),2));
figure; errorbar(tr(1,max_end-200:max_end+200)', avg_release, sem_release, '-k');
title(['Releases- avg +/- sem - n =' num2str(length(good_release_trials))])
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_release_avg.eps'], '-depsc');
print([dest '_release_avg.pdf'], '-dpdf');
savefig([dest '_release_avg.fig']);

first_third =1:round(ntrials/3);
last_third = round(2*ntrials/3):ntrials;
first_good_press_trials = intersect(first_third,good_press_trials);
last_good_press_trials = intersect(last_third,good_press_trials);
first_good_release_trials = intersect(first_third,good_release_trials);
last_good_release_trials = intersect(last_third,good_release_trials);

avg_first_third_press = nanmean(lever_press_mm(max_start-200:max_start+200,first_good_press_trials),2);
sem_first_third_press = nanstd(lever_press_mm(max_start-200:max_start+200,first_good_press_trials),[],2)./sqrt(sum(~isnan(lever_press_mm(max_start-200:max_start+200,first_good_press_trials)),2));
avg_last_third_press = nanmean(lever_press_mm(max_start-200:max_start+200,last_good_press_trials),2);
sem_last_third_press = nanstd(lever_press_mm(max_start-200:max_start+200,last_good_press_trials),[],2)./sqrt(sum(~isnan(lever_press_mm(max_start-200:max_start+200,last_good_press_trials)),2));
figure; errorbar(tp(1,max_start-200:max_start+200)', avg_first_third_press, sem_first_third_press, '-b');
hold on; errorbar(tp(1,max_start-200:max_start+200)', avg_last_third_press, sem_last_third_press, '-g');
title(['Avg Press- Blue: first thirdn = ' num2str(length(first_good_press_trials)) '; Green: last third n = ' num2str(length(last_good_press_trials))])
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_press_thirds.eps'], '-depsc');
print([dest '_press_thirds.pdf'], '-dpdf');
savefig([dest '_press_thirds.fig']);

avg_first_third_release = nanmean(lever_release_mm(max_end-200:max_end+200,first_good_release_trials),2);
sem_first_third_release = nanstd(lever_release_mm(max_end-200:max_end+200,first_good_release_trials),[],2)./sqrt(sum(~isnan(lever_release_mm(max_end-200:max_end+200,first_good_release_trials)),2));
avg_last_third_release = nanmean(lever_release_mm(max_end-200:max_end+200,last_good_release_trials),2);
sem_last_third_release = nanstd(lever_release_mm(max_end-200:max_end+200,last_good_release_trials),[],2)./sqrt(sum(~isnan(lever_release_mm(max_end-200:max_end+200,last_good_release_trials)),2));
figure; errorbar(tr(1,max_end-200:max_end+200)', avg_first_third_release, sem_first_third_release, '-b');
hold on; errorbar(tr(1,max_end-200:max_end+200)', avg_last_third_release, sem_last_third_release, '-g');
title(['Avg Release- Blue: first third n = ' num2str(length(first_good_release_trials)) '; Green: last third n = ' num2str(length(last_good_release_trials))])
xlabel('Time (ms)')
ylabel('Lever position (mm)')
print([dest '_release_thirds.eps'], '-depsc');
print([dest '_release_thirds.pdf'], '-dpdf');
savefig([dest '_release_thirds.fig']);

save([dest '_lever.mat'], 'good_press_trials', 'good_release_trials', 'lever_press_mm', 'lever_release_mm', 'max_start', 'max_end');