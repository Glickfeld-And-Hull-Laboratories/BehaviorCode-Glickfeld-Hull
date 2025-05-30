clear all
close all

%% Load behavior data
%mouse = "i2174";
%filenames = ["data-i2174-250225-1539.mat", "data-i2174-250226-1037.mat", "data-i2174-250228-1047.mat", "data-i2174-250301-1312.mat"]

%mouse = "i2182";
%filenames = ["data-i2182-250304-1456.mat", "data-i2182-250305-1055.mat", "data-i2182-250306-1140.mat", "data-i2182-250307-0944"]

%mouse = "i2179";
%filenames = ["data-i2179-250304-1459.mat", "data-i2179-250305-1059.mat", "data-i2179-250306-1144.mat", "data-i2179-250307-0948.mat"]';

%mouse = "i2140"
%filenames = ["data-i2140-240904-1152.mat", "data-i2140-240905-1151.mat", "data-i2140-240906-1055.mat", "data-i2140-240907-1130.mat"]

mouse = "i2167"
filenames = ["data-i2167-241211-1203.mat", "data-i2167-241213-1051.mat", "data-i2167-241214-1130.mat", "data-i2167-241215-1035.mat"]

dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
n_days = length(filenames);

for i = 1:length(filenames)
    load(fullfile(dir,filenames(i)));
    temp(i) = input;
end

%% Reaction time comparison (all trials)

react_time = [];
day = [];
dart = [];
background = [];

for i = 1:n_days
    i_react_time = double(cell2mat(temp(i).reactTimesMs));
    i_background = double(cell2mat(temp(i).tBlock2TrialNumber));

    n_trials = length(i_react_time);
    i_day = repelem(i, n_trials);
    i_dart = repelem(ceil(i/2), n_trials);

    react_time = [react_time, i_react_time];
    day = [day, i_day];
    dart = [dart, i_dart];
    background = [background, i_background];
end

for i = 1:length(dart)
    if dart(i) == 1
        dart_cat(i) = "pre-DART";
    elseif dart(i) == 2
        dart_cat(i) = "post-DART";
    end
end

for i = 1:length(background)
    if background(i) == 0
        background_cat(i) = "iso";
    elseif background(i) == 1
        background_cat(i) = "cross";
    end
end


fig1 = figure(1)
boxplot(react_time, {dart_cat, background_cat})
xlabel("Condition")
ylabel("Reaction time(ms)")
title([mouse, "Reaction times (all trials)"])

react_time_anova = anova({dart_cat', background_cat'}, react_time', FactorNames = ["DART", "Background"], ModelSpecification = "full")

%% Reaction time comparison (false alarms excluded)

outcome = [];

for i = 1:n_days
    i_outcome = string(temp(i).trialOutcomeCell);
    outcome = [outcome, i_outcome];
end

hits = strcmp(outcome,'success'); %Index of successes (ones)
misses = strcmp(outcome,'ignore'); %Index of misses (ones)
FAs = strcmp(outcome,'failure'); %Index of false alarms (ones)

hitIx = find(hits);
missIx = find(misses);
FAIx = find(FAs);

notFAs = hits + misses;

react_time_notFA =  react_time(find(notFAs));
background_cat_notFA = background_cat(find(notFAs));
day_notFA = day(find(notFAs));
dart_cat_notFA = dart_cat(find(notFAs));

fig2 = figure(2)
boxplot(react_time_notFA, {dart_cat_notFA, background_cat_notFA})
xlabel("Condition")
ylabel("Reaction time(ms)")
title([mouse, "Reaction times (false alarms excluded)"])

react_time_notFA_anova = anova({dart_cat_notFA', background_cat_notFA'}, react_time_notFA', FactorNames = ["DART", "Background"], ModelSpecification = "full")


%% Reaction time comparison (only hits)

react_time_hit = react_time(find(hits));
background_cat_hit = background_cat(find(hits));
day_hit = day(find(hits));
dart_cat_hit = dart_cat(find(hits));

fig3 = figure(3)
boxplot(react_time_hit, {dart_cat_hit, background_cat_hit})
axis square
xlabel("Condition")
ylabel("Reaction time(ms)")
title([mouse, "Reaction times (hits only)"])

react_time_hit_anova = anova({dart_cat_hit', background_cat_hit'}, react_time_hit', FactorNames = ["DART", "Background"], ModelSpecification = "full")

%{
figure(4)
histogram(react_time_hit)
%}


%% Psychometric curve (iso in one panel, cross in the other)

contrast_list = [];

for i = 1:n_days
    i_contrast_list = double(cell2mat(temp(i).gratingContrast));
    contrast_list = [contrast_list, i_contrast_list];
end


target_dart_condition = ["pre-DART", "post-DART", "pre-DART", "post-DART"];
target_background = ["iso", "iso", "cross", "cross"];


for i = 1:4
    x.condition = target_dart_condition(i) + "_" + target_background(i);
    x.trials = intersect(find(dart_cat == target_dart_condition(i)), find(background_cat == target_background(i)));
    x.valid_cons = contrast_list(x.trials);
    x.unique_cons = unique(x.valid_cons);
    x.pctcorr = zeros(length(x.unique_cons),1);
    x.pctFA = zeros(length(x.unique_cons),1);

    for j = 1:length(x.unique_cons)
        con = x.unique_cons(j);
        trial_list = x.trials(x.valid_cons == con);
        hit_list = intersect(trial_list, hitIx);
        miss_list = intersect(trial_list, missIx);
        FA_list = intersect(trial_list, FAIx);
        hit_num = length(hit_list);
        miss_num = length(miss_list);
        FA_num = length(FA_list);
        hit_rate = hit_num/(hit_num + miss_num);
        x.pctcorr(j) = hit_rate;
        FA_rate = FA_num/length(trial_list);
        x.pctFA(j) = FA_rate;
    end

    outcome_info(i) = x;

end
        
fig4 = figure(4)
sgtitle([mouse, "hit rate over contrast"])
subplot(1,2,1)
title("Iso")
axis square
hold on
plot(outcome_info(1).unique_cons, outcome_info(1).pctcorr, '-ok')
plot(outcome_info(2).unique_cons, outcome_info(2).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')
subplot(1,2,2)
title("Cross")
axis square
hold on
plot(outcome_info(3).unique_cons, outcome_info(3).pctcorr, '-ok')
plot(outcome_info(4).unique_cons, outcome_info(4).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')

%% Compare FA rates between days

FA_rates = zeros(1,n_days);

for i = 1:n_days
    n_total_trials = numel(find(day == i));
    n_FA = numel(intersect(find(day == i), find(FAs)));
    FA_rates(i) = n_FA/n_total_trials;
end


fig5 = figure(5)
bar(FA_rates)
xlabel('Experimental day')
ylabel('False alarm rate')
ylim([0,1])
title([mouse, "False alarm rates"])

%% FA rate by contrast

        
fig6 = figure(6)
sgtitle([mouse, "false alarm rate over contrast"])
subplot(1,2,1)
title("Iso")
axis square
hold on
plot(outcome_info(1).unique_cons, outcome_info(1).pctFA, '-ok')
plot(outcome_info(2).unique_cons, outcome_info(2).pctFA, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('False alarm rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')
subplot(1,2,2)
title("Cross")
axis square
hold on
plot(outcome_info(3).unique_cons, outcome_info(3).pctFA, '-ok')
plot(outcome_info(4).unique_cons, outcome_info(4).pctFA, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('False alarm rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')

%% Compare overall number of trials

trial_nums = zeros(1,n_days);

for i = 1:n_days
    n_total_trials = numel(find(day == i));
    trial_nums(i) = n_total_trials;
end

fig7 = figure(7)
bar(trial_nums)
xlabel('Experimental day')
ylabel('Total number of trials')
ylim([0,600])
title([mouse, "Total trial numbers"])

%% History dependence

%% Fit cumulative weibull to psychommetric curves (assisted by chatGPT)

ft = fittype('1 - exp(-(x / lambda)^beta)', 'independent', 'x', 'dependent', 'y');

for i = 1:length(outcome_info)
    opts = fitoptions(ft);
    opts.StartPoint = [max(outcome_info(i).unique_cons), 1];  % Initial guess: [lambda, beta]
    opts.Upper = [max(outcome_info(i).unique_cons) * 10, 1000];
    opts.Lower = [0, 0.0001];
    [outcome_info(i).weibullFit, outcome_info(i).gof] = fit(outcome_info(i).unique_cons(:), outcome_info(i).pctcorr(:), ft, opts);
end

fig8 = figure(8);
sgtitle([mouse, "hit rate over contrast"])
subplot(1,2,1)
title("Iso")
axis square
hold on
plot(outcome_info(1).unique_cons, outcome_info(1).pctcorr, '-ok')
plot(outcome_info(2).unique_cons, outcome_info(2).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
xRange_1 = linspace(min(outcome_info(1).unique_cons), max(outcome_info(1).unique_cons), 100);
yRange_1 = 1 - exp(-(xRange_1 / outcome_info(1).weibullFit.lambda).^outcome_info(1).weibullFit.beta);;
plot(xRange_1, yRange_1, 'k:', 'LineWidth', 1);
xRange_2 = linspace(min(outcome_info(2).unique_cons), max(outcome_info(2).unique_cons), 100);
yRange_2 = 1 - exp(-(xRange_2 / outcome_info(2).weibullFit.lambda).^outcome_info(2).weibullFit.beta);;
plot(xRange_2, yRange_2, 'b:', 'LineWidth', 1);
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')
subplot(1,2,2)
title("Cross")
axis square
hold on
plot(outcome_info(3).unique_cons, outcome_info(3).pctcorr, '-ok')
plot(outcome_info(4).unique_cons, outcome_info(4).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
xRange_3 = linspace(min(outcome_info(3).unique_cons), max(outcome_info(3).unique_cons), 100);
yRange_3 = 1 - exp(-(xRange_3 / outcome_info(3).weibullFit.lambda).^outcome_info(3).weibullFit.beta);;
plot(xRange_3, yRange_3, 'k:', 'LineWidth', 1);
xRange_4 = linspace(min(outcome_info(4).unique_cons), max(outcome_info(4).unique_cons), 100);
yRange_4 = 1 - exp(-(xRange_4 / outcome_info(4).weibullFit.lambda).^outcome_info(4).weibullFit.beta);;
plot(xRange_4, yRange_4, 'b:', 'LineWidth', 1);
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')


%% Psychometric curve (pre-dart in one panel, post-dart in the other) -- COMMENT OUT
%{
contrast_list = [];

for i = 1:n_days
    i_contrast_list = double(cell2mat(temp(i).gratingContrast));
    contrast_list = [contrast_list, i_contrast_list];
end


target_dart_condition = ["pre-DART", "pre-DART", "post-DART", "post-DART"];
target_background = ["iso", "cross", "iso", "cross"];


for i = 1:4
    x.condition = target_dart_condition(i) + "_" + target_background(i);
    x.trials = intersect(find(dart_cat == target_dart_condition(i)), find(background_cat == target_background(i)));
    x.valid_cons = contrast_list(x.trials);
    x.unique_cons = unique(x.valid_cons);
    x.pctcorr = zeros(length(x.unique_cons),1);
    x.pctFA = zeros(length(x.unique_cons),1);

    for j = 1:length(x.unique_cons)
        con = x.unique_cons(j);
        trial_list = x.trials(x.valid_cons == con);
        hit_list = intersect(trial_list, hitIx);
        miss_list = intersect(trial_list, missIx);
        FA_list = intersect(trial_list, FAIx);
        hit_num = length(hit_list);
        miss_num = length(miss_list);
        FA_num = length(FA_list);
        hit_rate = hit_num/(hit_num + miss_num);
        x.pctcorr(j) = hit_rate;
        FA_rate = FA_num/length(trial_list);
        x.pctFA(j) = FA_rate;
    end

    outcome_info_panels_flipped(i) = x;

end
        
fig90 = figure(90)
sgtitle([mouse, "hit rate over contrast"])
subplot(1,2,1)
title("pre-DART")
axis square
hold on
plot(outcome_info_panels_flipped(1).unique_cons, outcome_info_panels_flipped(1).pctcorr, '-ok')
plot(outcome_info_panels_flipped(2).unique_cons, outcome_info_panels_flipped(2).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
legend('Iso', 'Cross', 'Location', 'southeast')
subplot(1,2,2)
title("post-DART")
axis square
hold on
plot(outcome_info_panels_flipped(3).unique_cons, outcome_info_panels_flipped(3).pctcorr, '-ok')
plot(outcome_info_panels_flipped(4).unique_cons, outcome_info_panels_flipped(4).pctcorr, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('Hit rate')
legend('Iso', 'Cross', 'Location', 'southeast')

fig91 = figure(91)
sgtitle([mouse, "false alarm rate over contrast"])
subplot(1,2,1)
title("Iso")
axis square
hold on
plot(outcome_info_panels_flipped(1).unique_cons, outcome_info_panels_flipped(1).pctFA, '-ok')
plot(outcome_info_panels_flipped(2).unique_cons, outcome_info_panels_flipped(2).pctFA, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('False alarm rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')
subplot(1,2,2)
title("Cross")
axis square
hold on
plot(outcome_info_panels_flipped(3).unique_cons, outcome_info_panels_flipped(3).pctFA, '-ok')
plot(outcome_info_panels_flipped(4).unique_cons, outcome_info_panels_flipped(4).pctFA, '-ob')
ylim([0,1])
xscale log
xlabel('Contrast')
ylabel('False alarm rate')
legend('Pre-DART', 'Post-DART', 'Location', 'southeast')

%}


%% Export
%{
savedir = fullfile('Z:\All_Staff\home\ACh\Analysis\Behavior', mouse)

saveas(fig3, fullfile(savedir, 'reaction_times.jpg'))
saveas(fig4, fullfile(savedir, 'psychometric_curves.jpg'))
saveas(fig5, fullfile(savedir, 'FA_rates_between_days.jpg'))
saveas(fig6, fullfile(savedir, 'FA_rates_by_condition.jpg'))
saveas(fig7, fullfile(savedir, 'total_trial_numbers.jpg'))
%saveas(fig90, fullfile(savedir, 'psychometric_curves_panels_flipped.jpg'))
%saveas(fig91, fullfile(savedir, 'FA_rates_by_condition_panels_flipped.jpg'))

%}



    










