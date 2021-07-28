% This script sorts the data from the test2AFCDatasetGenerator.m task and generates a GLM for the
% control trials
clear all

%% params
t_spo = 2;
t_max = 1;
d_max = 100;
d_spo = 0.37;
n = 1:4;

%calculate all targets and differences
t =  fliplr(t_max./(2.^((n-1)./(t_spo))));
d =  fliplr(d_max./(2.^((n-1)./(d_spo)))+1);

%random trial generator
ntrials = 2516;
x_t = randi(n(end),[1 ntrials]); %which target [1:4]
x_d = randi(n(end),[1 ntrials]); %which diff [1:4]
x_s = randi([0 1],[1 ntrials]); %target location: 0-left 1-right

x_tcon = t(x_t); %actual target con
x_dcon = x_tcon./d(x_d); %actual distractor con
r_con = x_tcon;
r_con(find(x_s==0)) = x_dcon(find(x_s==0)); %actual right con
l_con = x_tcon;
l_con(find(x_s==1)) = x_dcon(find(x_s==1)); %actual left con

%% contrast discrim - good mouse we want
%probability of correct choice by contrast diff
%d_p = [0.6 0.7 0.85 1]; arbitrary
d_p = [0.4777 0.8024 0.9733 0.9711]; % based on i459's actual prob right
y_s = zeros(1,ntrials);
for i = 1:ntrials
    if rand(1)<=d_p(x_d(i))
        y_s(i) = x_s(i);
    else
        y_s(i) = abs(x_s(i)-1);
    end
end

%% Randomize and separate data

nTrialsB1 = ntrials;


trialType = x_s';

% trialType -1 is a left trial, 1 is a right trial 
trialType(trialType == 0) = -1;

tGratingContrast = x_tcon';
dGratingContrast = x_dcon';
tRightResponse = y_s';


all_variables_B1 = [trialType, tGratingContrast, dGratingContrast, tRightResponse];

s = rng('default');

rando = randperm(nTrialsB1); %arrary of random positions
rand_dataB1 = zeros(size(all_variables_B1));
for k = 1:nTrialsB1
    rand_dataB1(k,:) = all_variables_B1(rando(k),:);
end


train_dataB1 = rand_dataB1(1:floor(nTrialsB1/2),:);
test_dataB1 = rand_dataB1(1+floor(nTrialsB1/2):nTrialsB1,:);

%% Training

% make table with binary labels for each parameter

trialType = train_dataB1(:,1);

contrast_ratios = round(abs(train_dataB1(:,2)./train_dataB1(:,3)),3); % target/distractor grating contrast 

contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);

% % Ratio processing
% if length(unique(contrast_ratios)) > 8  
%     contrast_ratios = wMeanRatios(contrast_ratios);
%     unique_contrast_ratios = unique(contrast_ratios);
% end

unique_contrast_ratios = unique(contrast_ratios);
zRatios = zscore(contrast_ratios);
unique_zRatios = unique(zRatios);

target_contrasts = round(abs(train_dataB1(:,2)),3);
target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
unique_target_contrasts = unique(target_contrasts);
zContrasts = zscore(target_contrasts);
unique_zContrasts = unique(zContrasts);

RatiosContrasts = contrast_ratios.*target_contrasts;
zRatiosContrasts = zscore(RatiosContrasts);
unique_zRatiosContrasts = unique(zRatiosContrasts);


mean_zScores = mean([zRatios zContrasts zRatiosContrasts]); % should be close to 0


tRightResponse = train_dataB1(:,4); 

parameter_table = table(tRightResponse, zRatios, zContrasts,  zRatiosContrasts);

%% compute the model fit 

model_spec_2AFC = 'tRightResponse ~ zRatios + zContrasts + zRatiosContrasts';

model_2AFC = fitglm(parameter_table, model_spec_2AFC, 'distribution', 'binomial')

clone_B1nTrials = length(tRightResponse);

clone_model_2AFC = model_2AFC;

save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
        'clone_modelInfo.mat'],'clone_model_2AFC', 'clone_B1nTrials');

% figure
% plotResiduals(model_2AFC);

%% get block 1 zscored weights (will be used in later figure with b2 weights)

ratio_weights = model_2AFC.Coefficients.Estimate(2);
ratio_weights_se = model_2AFC.Coefficients.SE(2);
contrast_weights = model_2AFC.Coefficients.Estimate(3);
contrast_weights_se = model_2AFC.Coefficients.SE(3);
RatioContrast_weights = model_2AFC.Coefficients.Estimate(4);
RatioContrast_se = model_2AFC.Coefficients.SE(4);
bias = model_2AFC.Coefficients.Estimate(1);
bias_se = model_2AFC.Coefficients.SE(1);


r1 = ratio_weights;
r1_se = ratio_weights_se;
c1 = contrast_weights;
c1_se = contrast_weights_se;
b1 = bias;
b1_se = bias_se;
rtc1 = RatioContrast_weights;
rtc1_se = RatioContrast_se;

data = [r1 ; c1 ; rtc1 ; b1];
err = [r1_se ; c1_se ; rtc1_se ; b1_se ];
figure(2)
% categories = categorical({'Contrast Ratio', 'Target Contrast', 'Ratios*Contrasts', 'Bias'});
weights = [ratio_weights contrast_weights RatioContrast_weights bias];
deviations = [ratio_weights_se contrast_weights_se RatioContrast_se bias_se];
bar(data, 0.25)
% bar(categories, weights, 0.25);
box('off')
set(gca,'FontSize',10)
set(gca, 'XTickLabel', {'Contrast Ratio' 'Target Contrast' 'Ratio*Contrast' 'Bias'});
ylabel('weight')
hold on

err1 = errorbar(data, err, 'linestyle', 'none', 'color', 'k');


% err1 = errorbar(categories, weights, deviations, 'linestyle', 'none', 'color', 'k');
er.color = [0 0 0];                            
er.lineStyle = 'none'; 

% title(['Fake Mouse %Right Model'])
title(['Clone %Right Model'])


hold off

saveas(gcf, 'clone_weights.png')

% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% 'clone_weights.pdf'],'-dpdf','-bestfit');


%% Testing

% get the contrast ratio of each trial and unique contrasts
contrast_ratios_nogo = round(abs(test_dataB1(:,2)./test_dataB1(:,3)),3); % target/distractor grating contrast 

% "Nogo" trials have equal contrast ratios. We do not want to include them,
% so remove them. (this trial type is only found in  in early mice)
test_dataB1 = test_dataB1(contrast_ratios_nogo ~= 1,:); 

% make table with binary labels for each parameter

trialType = test_dataB1(:,1);

contrast_ratios = round(abs(test_dataB1(:,2)./test_dataB1(:,3)),3); % target/distractor grating contrast 

contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);

% % Ratio processing
% if length(unique(contrast_ratios)) > 8  
%     contrast_ratios = wMeanRatios(contrast_ratios);
%     unique_contrast_ratios = unique(contrast_ratios);
% end

unique_contrast_ratios = unique(contrast_ratios);
zRatios = zscore(contrast_ratios);
unique_zRatios = unique(zRatios);

target_contrasts = round(abs(test_dataB1(:,2)),3);
target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
unique_target_contrasts = unique(target_contrasts);
zContrasts = zscore(target_contrasts);
unique_zContrasts = unique(zContrasts);

RatiosContrasts = contrast_ratios.*target_contrasts;
zRatiosContrasts = zscore(RatiosContrasts);
unique_zRatiosContrasts = unique(zRatiosContrasts);

mean_zScores = mean([zRatios zContrasts zRatiosContrasts]);


tRightResponse = test_dataB1(:,4); 

test_parameter_table = table(tRightResponse, zRatios, zContrasts, zRatiosContrasts);


%% Use predict function and trained model to predict responses for testing data

[predict_tRightResponse, predictCI] = predict(model_2AFC, test_parameter_table);

%% Set threshold for responses

%use instead of hard threshold

nTrials = length(test_dataB1);

out = predict_tRightResponse;
right_thresh = zeros(1,nTrials);
labels = zeros(1,nTrials);
for i = 1:nTrials
  labels(i) = rand(1);
  if labels(i)<out(i)
    right_thresh(i) = 1;
  end
end

right_thresh = right_thresh';
% predict_parameter_table = table(right_thresh, tRatios, tContrasts);

 
 % compare test and predict (% right)

nActualRightResponseB1 = [];
totalTrials = [];
for i = 1:length(unique_contrast_ratios)
    nActualRightResponseB1(i) = sum(tRightResponse(contrast_ratios == unique_contrast_ratios(i))); %correct right responses per stimulus level
    totalTrials(i) = sum(contrast_ratios == unique_contrast_ratios(i)); %total right trials per stimulus level
end

nPredictRightResponseB1 = [];

for i = 1:length(unique_contrast_ratios)
    nPredictRightResponseB1(i) = sum(right_thresh(contrast_ratios == unique_contrast_ratios(i))); % predicted right responses per stimulus level
end


 %%       

%find success per contrast ratio and totals for binomial fit

[phat, pci] = binofit(nActualRightResponseB1, totalTrials);

neg = [];
pos =[];

for i= 1:length(phat)
    neg(i) = phat(i) - pci(i,1);
    pos(i) = pci(i,2) - phat(i);
end


%% make plot for test data...(actual)

xValues(1:4) = -1./unique_contrast_ratios(1:4);
xValues(5:8) = unique_contrast_ratios(5:8);

figure(1)
B1TestPlot = errorbar(xValues, phat, neg, pos, 'Linewidth', 1.5, 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', 'Test');
hold on 


% vH = plot([0,0],[0,1]);
% set(vH, 'Color', 'g');
% % set(gca,'xscale','log')
% xlabel('Contrast Ratio (R/L)')
% ylabel('% Right')
% grid on
% xlim([-105 105])
% title(['Control Trials: n = '  num2str(B1nTrials)])
% 
% hold off


%%
%...and with predicted data (model)
 
[mphat, mpci] = binofit(nPredictRightResponseB1, totalTrials);

neg = [];
pos =[];

for i= 1:length(mphat)
    neg(i) = mphat(i) - mpci(i,1);
    pos(i) = mpci(i,2) - mphat(i);
end


YB1 = nActualRightResponseB1;
UB1 = nPredictRightResponseB1;

EVB = 1-(var(YB1-UB1)/var(YB1))



figure(1)
B1modelPlot = errorbar(xValues, mphat, neg, pos, 'Linewidth', 1.5, 'Marker', 's', 'MarkerSize', 8, 'Color', 'k', 'DisplayName', 'Testing' );

hold on
vH = plot([0,0],[0,1]);
set(vH, 'Color', 'g');
xlabel('Contrast Ratio (R/L)')
set(gca,'Xscale','log')
ylabel('% Right')
grid on
xlim([0 105])
ylim([0 1])
title(['Clone, Control: n = '  num2str(clone_B1nTrials) '; EV = ' num2str(EVB)])

hold off

legend([B1TestPlot B1modelPlot], {'Actual', 'Model'}, 'Location', 'northwest')

saveas(gcf, 'clone_psych.png')