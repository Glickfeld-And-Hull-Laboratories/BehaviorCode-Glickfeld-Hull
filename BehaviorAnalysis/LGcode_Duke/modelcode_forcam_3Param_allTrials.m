% This script sorts the data from the 2AFC task and generates a GLM for the
% control trials

%Original author: Robin Blazing; edits Camaron Mangham
%%%%%%%%%%%%%%

%%
clear all

%% Load input data (ex: "i459goodData")
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim';
mouse = 'i459';
load(fullfile(data_path,[mouse 'goodData.mat']));


%% Process training data into parameter tables (design matrix) for use in GLM

% concatenate all entries in structure fields to single arrays 
SIx = [separated_input_data.correctIx]'; %successes
FIx = [separated_input_data.incorrectIx]'; %failures
tLeftTrial = [separated_input_data.tLeftTrial]'; % trial type (right or left)
tGratingContrast = [separated_input_data.tGratingContrast]'; % target grating contrast
dGratingContrast = [separated_input_data.dGratingContrast]'; % distractor grating contrast
tRightResponse = [separated_input_data.tRightResponse]'; % response type (right or left)
tIgnore = [separated_input_data.ignoreIx]'; %index of ignored trials
tBlock2 = [separated_input_data.tBlock2TrialNumber]'; %index of optogenetic stimulation trials
rightGratingContrast = [separated_input_data.rightGratingContrast]';
leftGratingContrast = [separated_input_data.leftGratingContrast]';


% make sure object class is double
SIx = double(SIx);
FIx = double(FIx);
tLeftTrial = double(celleqel2mat_padded(tLeftTrial));
tGratingContrast = double(celleqel2mat_padded(tGratingContrast));
dGratingContrast = double(celleqel2mat_padded(dGratingContrast));
tRightResponse = double(celleqel2mat_padded(tRightResponse));
tIgnore = double(tIgnore);
tBlock2 = double(celleqel2mat_padded(tBlock2));
rightGratingContrast = double(celleqel2mat_padded(rightGratingContrast));
leftGratingContrast = double(celleqel2mat_padded(leftGratingContrast));

% get number of trials
nTrials = length(tLeftTrial); 

% % make an array for whether the trial was a left choice (save for later?)
tLeftChoice = zeros(1,nTrials); 
tLeftChoice(SIx == 1 & tLeftTrial == 1) = 1;
tLeftChoice(FIx == 1 & tLeftTrial == 0) = 1;

% new variable: trialType -1 is a left trial, 1 is a right trial 
trialType(tLeftTrial == 1) = -1;
trialType(tLeftTrial == 0) = 1;


 sizeTT = size(trialType);
    if sizeTT(1) == 1
        trialType = trialType';
    end


%concatenate variables into array so that they are easier to manipulate 
all_variables = [SIx, FIx, trialType, tGratingContrast, dGratingContrast, tRightResponse, tIgnore, tBlock2,...
    rightGratingContrast, leftGratingContrast]; 

% remove ignores and block 2
all_variables_B1 = all_variables(tIgnore == 0 & tBlock2 == 0, :);
all_variables_B2 = all_variables(tIgnore == 0 & tBlock2 == 1, :);

%% Randomize and separate data

nTrialsB1 = length(all_variables_B1);

s = rng('default');

rando = randperm(nTrialsB1); %arrary of random positions
rand_dataB1 = zeros(size(all_variables_B1));
for k = 1:nTrialsB1
    rand_dataB1(k,:) = all_variables_B1(rando(k),:);
end


train_dataB1 = rand_dataB1(1:floor(nTrialsB1/2),:);
test_dataB1 = rand_dataB1(1+floor(nTrialsB1/2):nTrialsB1,:);

%Block 2

nTrialsB2 = length(all_variables_B2);

s = rng('default');

rando = randperm(nTrialsB2); %arrary of random positions
rand_dataB2 = zeros(size(all_variables_B2));
for k = 1:nTrialsB2
    rand_dataB2(k,:) = all_variables_B2(rando(k),:);
end


train_dataB2 = rand_dataB2(1:floor(nTrialsB2/2),:);
test_dataB2 = rand_dataB2(1+floor(nTrialsB2/2):nTrialsB2,:);

%% Training

% get the contrast ratio of each trial and unique contrasts
contrast_ratios_nogo = round(abs(train_dataB1(:,4)./train_dataB1(:,5)),3); % target/distractor grating contrast 

% "Nogo" trials have equal contrast ratios. We do not want to include them,
% so remove them. (this trial type is only found in  in early mice)
train_dataB1 = train_dataB1(contrast_ratios_nogo ~= 1,:); 

% make table with binary labels for each parameter

trialType = train_dataB1(:,3);

contrast_ratios = round(abs(train_dataB1(:,4)./train_dataB1(:,5)),3); % target/distractor grating contrast 

contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);

% Ratio processing
if length(unique(contrast_ratios)) > 8  
    contrast_ratios = wMeanRatios(contrast_ratios);
    unique_contrast_ratios = unique(contrast_ratios);
end

unique_contrast_ratios = unique(contrast_ratios);
tRatios = zscore(contrast_ratios);
unique_tRatios = unique(tRatios);

target_contrasts = round(abs(train_dataB1(:,4)),3);
target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
unique_target_contrasts = unique(target_contrasts);
tContrasts = zscore(target_contrasts);
unique_tContrast = unique(tContrasts);

tRatios_t_tContrasts = contrast_ratios .* target_contrasts;
tRatios_t_tContrasts = zscore(tRatios_t_tContrasts);

mean_zScores = mean([tRatios tContrasts tRatios_t_tContrasts]);


tRightResponse = train_dataB1(:,6); 

parameter_table = table(tRightResponse, tRatios, tContrasts, tRatios_t_tContrasts);


%% Set up parameter tables
% 
% % tRightResponse = training_data(:,1);
% % contrast_ratios = training_data(:,3);
% % unique_contrast_ratios = unique(contrast_ratios);
% 
% parameter_table = table(training_data(:,1), training_data(:,2), training_data(:,3), 'VariableNames',...
% {'tRightResponse', 'tRatios', 'tContrasts'});
% 
% test_parameter_table = table(test_data(:,1), test_data(:,2), test_data(:,3), 'VariableNames',...
% {'tRightResponse', 'tRatios', 'tContrasts'});


%% compute the model fit 

model_spec_2AFC = 'tRightResponse ~ tRatios + tContrasts + tRatios_t_tContrasts';

model_2AFC = fitglm(parameter_table, model_spec_2AFC, 'distribution', 'binomial')

B1nTrials = length(tRightResponse);

save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
        mouse '_modelInfo.mat'],'model_2AFC', 'B1nTrials');

% figure
% plotResiduals(model_2AFC);

%% get block 1 zscored weights (will be used in later figure with b2 weights)

ratio_weights = model_2AFC.Coefficients.Estimate(2);
ratio_weights_se = model_2AFC.Coefficients.SE(2);
contrast_weights = model_2AFC.Coefficients.Estimate(3);
contrast_weights_se = model_2AFC.Coefficients.SE(3);
bias = model_2AFC.Coefficients.Estimate(1);
bias_se = model_2AFC.Coefficients.SE(1);
rtc_weights = model_2AFC.Coefficients.Estimate(4);
rtc_weights_se = model_2AFC.Coefficients.SE(4);


r1 = ratio_weights;
r1_se = ratio_weights_se;
c1 = contrast_weights;
c1_se = contrast_weights_se;
b1 = bias;
b1_se = bias_se;
rtc1 = rtc_weights;
rtc1_se = rtc_weights_se;


data = [r1 ; c1 ; rtc1 ; b1];
err = [r1_se ; c1_se ; rtc1_se ; b1_se ];
figure(2)
% categories = categorical({'Contrast Ratio', 'Target Contrast', 'Ratios*Contrasts', 'Bias'});
% weights = [ratio_weights contrast_weights RatioContrast_weights bias];
% deviations = [ratio_weights_se contrast_weights_se RatioContrast_se bias_se];
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
title([mouse ' %Right Model'])

saveas(gcf, [mouse '_weights.png'])
% 
% 
% hold off

% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% mouse '_weights.pdf'],'-dpdf','-bestfit');


 %% plot training data   

% %find sucssessful right responses per contrast ratio and totals for binomial fit
% 
% nRightResponse = [];
% totalTrials = [];
% for i = 1:length(unique_contrast_ratios)
%     nRightResponse(i) = sum(tRightResponse(contrast_ratios == unique_contrast_ratios(i)));
%     totalTrials(i) = sum((contrast_ratios == unique_contrast_ratios(i)));
% end
% 
% [trainphat, trainpci] = binofit(nRightResponse, totalTrials);
% 
% neg = [];
% pos =[];
% 
% for i= 1:length(trainphat)
%     neg(i) = trainphat(i) - trainpci(i,1);
%     pos(i) = trainpci(i,2) - trainphat(i);
% end
% 
% figure(1)
% subplot(1,2,1);
% TrainingPlot = errorbar(unique_contrast_ratios, trainphat, neg, pos, 'Linewidth', 1.5, 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', 'Training');
% hold on 
% 
% 
% vH = plot([0,0],[0,1]);
% set(vH, 'Color', 'g');
% set(gca,'xscale','log')
% 
% xlabel('Contrast Ratio (R/L)')
% ylabel('% Right')
% grid on
% title('% Right by Contrast')
% legend([TrainingPlot], {'Training'}, 'Location', 'northwest')

%% Testing

% get the contrast ratio of each trial and unique contrasts
contrast_ratios_nogo = round(abs(test_dataB1(:,4)./test_dataB1(:,5)),3); % target/distractor grating contrast 

% "Nogo" trials have equal contrast ratios. We do not want to include them,
% so remove them. (this trial type is only found in  in early mice)
test_dataB1 = test_dataB1(contrast_ratios_nogo ~= 1,:); 

% make table with binary labels for each parameter

trialType = test_dataB1(:,3);

contrast_ratios = round(abs(test_dataB1(:,4)./test_dataB1(:,5)),3); % target/distractor grating contrast 

contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);

% Ratio processing
if length(unique(contrast_ratios)) > 8  
    contrast_ratios = wMeanRatios(contrast_ratios);
    unique_contrast_ratios = unique(contrast_ratios);
end

unique_contrast_ratios = unique(contrast_ratios);
tRatios = zscore(contrast_ratios);
unique_tRatios = unique(tRatios);

target_contrasts = round(abs(test_dataB1(:,4)),3);
target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
unique_target_contrasts = unique(target_contrasts);
tContrasts = zscore(target_contrasts);
unique_tContrast = unique(tContrasts);

tRatios_t_tContrasts = contrast_ratios .* target_contrasts;
tRatios_t_tContrasts = zscore(tRatios_t_tContrasts);

mean_zScores = mean([tRatios tContrasts tRatios_t_tContrasts]);


tRightResponse = test_dataB1(:,6); 

test_parameter_table = table(tRightResponse, tRatios, tContrasts, tRatios_t_tContrasts);


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
figure(1)

xValues(1:4) = -1./unique_contrast_ratios(1:4);
xValues(5:8) = unique_contrast_ratios(5:8);

B1TestPlot = errorbar(xValues, phat, neg, pos, 'Linewidth', 1.5, 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', 'Test');
hold on 

% set(gca,'xscale','log')
% vH = plot([0,0],[0,1]);
% set(vH, 'Color', 'g');
% set(gca,'xscale','log')
% xlabel('Contrast Ratio (R/L)')
% ylabel('% Right')
% grid on
% title('% Right by Contrast')


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

EVB1 = 1-(var(YB1-UB1)/var(YB1))

figure(1)
B1modelPlot = errorbar(xValues, mphat, neg, pos, 'Linewidth', 1.5, 'Marker', 's', 'MarkerSize', 8, 'Color', 'k', 'DisplayName', 'Testing' );

hold on
vH = plot([0,0],[0,1]);
set(vH, 'Color', 'g');
set(gca,'xscale','log')
xlabel('Contrast Ratio (R/L)')
ylabel('% Right')
grid on
title([mouse '- Control: n = '  num2str(B1nTrials)  '; EV = ' num2str(EVB1)])
xlim([0 105])
ylim([0 1])
hold off

legend([B1TestPlot B1modelPlot], {'Actual', 'Model'}, 'Location', 'northwest')

saveas(gcf, [mouse '_psych.png'])

% vH = plot([0,0],[0,1]);
% set(vH, 'Color', 'g');
% set(gca,'xscale','log')
% 
% xlabel('Contrast Ratio (R/L)')
% ylabel('% Right')
% grid on
% title('% Right by Contrast')
% % legend([TestPlot modelPlot], {'Testing', 'Model'}, 'Location', 'northwest')

% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% mouse '_B1Model_summary.pdf'],'-dpdf','-bestfit');
% 
% saveas(gcf,[mouse '_B1Model_summary', '.png'])


%% START BLOCK 2
% 
% 
% 
% 
% %START BLOCK 2 >>>>>>>>>>>
% 
% %% Process training data into parameter tables (design matrix) for use in GLM
% 
% 
% % get the contrast ratio of each trial and unique contrasts
% contrast_ratios_nogo = round(abs(train_dataB2(:,4)./train_dataB2(:,5)),3); % target/distractor grating contrast 
% 
% % "Nogo" trials have equal contrast ratios. We do not want to include them,
% % so remove them. (this trial type is only found in  in early mice)
% train_dataB2 = train_dataB2(contrast_ratios_nogo ~= 1,:); 
% 
% % make table with binary labels for each parameter
% 
% trialType = train_dataB2(:,3);
% 
% contrast_ratios = round(abs(train_dataB2(:,4)./train_dataB2(:,5)),3); % target/distractor grating contrast 
% 
% contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);
% 
% % Ratio processing
% if length(unique(contrast_ratios)) > 8  
%     contrast_ratios = wMeanRatios(contrast_ratios);
%     unique_contrast_ratios = unique(contrast_ratios);
% end
% 
% unique_contrast_ratios = unique(contrast_ratios);
% tRatios = zscore(contrast_ratios);
% unique_tRatios = unique(tRatios);
% 
% target_contrasts = round(abs(train_dataB2(:,4)),3);
% target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
% unique_target_contrasts = unique(target_contrasts);
% tContrasts = zscore(target_contrasts);
% unique_tContrast = unique(tContrasts);
% 
% tRatios_t_tContrasts = contrast_ratios .* target_contrasts;
% tRatios_t_tContrasts = zscore(tRatios_t_tContrasts);
% 
% 
% mean_zScores = mean([tRatios tContrasts tRatios_t_tContrasts]);
% 
% 
% tRightResponse = train_dataB2(:,6); 
% 
% B2parameter_table = table(tRightResponse, tRatios, tContrasts, tRatios_t_tContrasts);
% 
% 
% 
% 
% %% compute the model fit and plot residuals
% 
% B2model_spec_2AFC = ['tRightResponse ~ tRatios + tContrasts + tRatios_t_tContrasts'];
% 
% B2model_2AFC = fitglm(B2parameter_table, B2model_spec_2AFC, 'distribution', 'binomial')
% 
% B2nTrials = length(tRightResponse);
% 
% % save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% %         mouse '_modelInfo.mat'], 'model_2AFC_tRatios', 'model_2AFC_tContrasts',...
% %         'model_2AFC', 'B2model_2AFC_tRatios', 'B2model_2AFC_tContrasts',...
% %         'B2model_2AFC', 'B1nTrials', 'B2nTrials');
% 
% % figure
% % plotResiduals(model_2AFC);
% 
% %% plot zscored weights (how about with multiple mice?)
% 
% ratio_weights = B2model_2AFC.Coefficients.Estimate(2);
% ratio_weights_se = B2model_2AFC.Coefficients.SE(2);
% contrast_weights = B2model_2AFC.Coefficients.Estimate(3);
% contrast_weights_se = B2model_2AFC.Coefficients.SE(3);
% bias = B2model_2AFC.Coefficients.Estimate(1);
% bias_se = B2model_2AFC.Coefficients.SE(1);
% rtc_weights = B2model_2AFC.Coefficients.Estimate(4);
% rtc_weights_se = B2model_2AFC.Coefficients.SE(4);
% 
% 
% 
% % figure(1)
% % hold on
% % categories = categorical({'Target Ratio', 'Target Contrast'});
% % categories = reordercats(categories, {});
% % weights = [ratio_weights contrast_weights];
% % deviations = [ratio_weights_se contrast_weights_se];
% % bar2 = bar(categories, weights, 0.25);
% % box('off')
% % set(gca,'FontSize',15)
% % ylabel('weight')
% % hold on
% % err2 = errorbar(categories, weights, deviations, deviations, 'linestyle', 'none', 'color', 'k');
% % er.color = [0 0 0];                            
% % er.lineStyle = 'none'; 
% % title([mouse ' Model'])
% % hold off
% 
% %Control and B2 on same bar graph
% 
% r2 = ratio_weights;
% r2_se = ratio_weights_se;
% c2 = contrast_weights;
% c2_se = contrast_weights_se;
% b2 = bias;
% b2_se = bias_se;
% rtc2 = rtc_weights;
% rtc2_se = rtc_weights_se;
% 
% data = [r1 r2; c1 c2; rtc1 rtc2; b1 b2];
% err = [r1_se r2_se; c1_se c2_se; rtc1_se rtc2_se; b1_se b2_se];
% figure(3)
% bar(data)
% set(gca, 'XTickLabel', {'Contrast Ratio' 'Target Contrast' 'Ratio*Contrast' 'Bias'});
% box('off')
% set(gca,'FontSize',10)
% hold on
% 
% ngroups = size(data,1);
% nbars = size(data, 2);
% groupwidth = min(0.8, nbars/(nbars +1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, data(:,i), err(:,i), 'linestyle', 'none', 'color', 'k');
% end
% 
% ylabel('weight')
% legend('Control', 'Led Stim');
% title([mouse ' % Right Model Weights'])
% hold off
% 
% % print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% % mouse '_weights_summary.pdf'],'-dpdf','-bestfit');
% % 
% % saveas(gcf,[mouse '_weights_summary', '.png'])
% 
% 
% %% Is the difference in weights significant (Paired t-test)? 
% 
% % ratio_weights = model_2AFC.Coefficients.Estimate(2);
% % B2ratio_weights = B2model_2AFC.Coefficients.Estimate(2);
% % 
% % contrast_weights = model_2AFC.Coefficients.Estimate(3);
% % B2contrast_weights = B2model_2AFC.Coefficients.Estimate(3);
% % 
% % % [ratio_h, ratio_p] = ttest(ratio_weights, B2ratio_weights)
% % % [contrast_h, contrast_p] = ttest(contrast_weights, B2contrast_weights)
% 
% 
% 
% %% find sucssessful right responses per contrast ratio and totals for binomial fit
% 
% nRightResponse = [];
% totalTrials = [];
% for i = 1:length(unique_contrast_ratios)
%     nRightResponse(i) = sum(tRightResponse(contrast_ratios == unique_contrast_ratios(i)));
%     totalTrials(i) = sum((contrast_ratios == unique_contrast_ratios(i)));
% end
% 
% [trainphat, trainpci] = binofit(nRightResponse, totalTrials);
% 
% neg = [];
% pos =[];
% 
% for i= 1:length(trainphat)
%     neg(i) = trainphat(i) - trainpci(i,1);
%     pos(i) = trainpci(i,2) - trainphat(i);
% end
% 
%  %% plot block2 training data (plot on B1 training data...) (uncessary?)
% 
% 
% % subplot(1,2,1);
% % TrainingPlot = errorbar(unique_contrast_ratios, trainphat, neg, pos, 'Linewidth', 1.5, 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', 'Training');
% % hold on 
% 
% 
% % vH = plot([0,0],[0,1]);
% % set(vH, 'Color', 'g');
% % set(gca,'xscale','log')
% % 
% % xlabel('Contrast Ratio (R/L)')
% % ylabel('% Right')
% % grid on
% % title('Stim % Right by Contrast')
% % legend([TrainingPlot], {'Training'}, 'Location', 'northwest')
% 
% 
% %% Process block 2 testing data into parameter tables (design matrix) for use in GLM
% 
% % get the contrast ratio of each trial and unique contrasts
% contrast_ratios_nogo = round(abs(test_dataB2(:,4)./test_dataB2(:,5)),3); % target/distractor grating contrast 
% 
% % "Nogo" trials have equal contrast ratios. We do not want to include them,
% % so remove them. (this trial type is only found in  in early mice)
% test_dataB2 = test_dataB2(contrast_ratios_nogo ~= 1,:); 
% 
% % make table with binary labels for each parameter
% 
% trialType = test_dataB2(:,3);
% 
% contrast_ratios = round(abs(test_dataB2(:,4)./test_dataB2(:,5)),3); % target/distractor grating contrast 
% 
% contrast_ratios(trialType == -1) = -1.*contrast_ratios(trialType == -1);
% 
% % Ratio processing
% if length(unique(contrast_ratios)) > 8  
%     contrast_ratios = wMeanRatios(contrast_ratios);
%     unique_contrast_ratios = unique(contrast_ratios);
% end
% 
% unique_contrast_ratios = unique(contrast_ratios);
% tRatios = zscore(contrast_ratios);
% unique_tRatios = unique(tRatios);
% 
% target_contrasts = round(abs(test_dataB2(:,4)),3);
% target_contrasts(trialType == -1) = -1.*target_contrasts(trialType == -1);
% unique_target_contrasts = unique(target_contrasts);
% tContrasts = zscore(target_contrasts);
% unique_tContrast = unique(tContrasts);
% 
% tRatios_t_tContrasts = contrast_ratios .* target_contrasts;
% tRatios_t_tContrasts = zscore(tRatios_t_tContrasts);
% 
% mean_zScores = mean([tRatios tContrasts tRatios_t_tContrasts]);
% 
% 
% tRightResponse = test_dataB2(:,6);  
% 
% test_parameter_table_tBlock2_alone = table(tRightResponse);
% test_parameter_table_tBlock2 = table(tRightResponse, tRatios, tContrasts, tRatios_t_tContrasts);
% 
% 
% %% Use predict function and trained model to predict responses for testing data
% 
% [predict_tRightResponse, predictCI] = predict(B2model_2AFC, test_parameter_table_tBlock2);
% 
% %% Set threshold for responses
% 
% %use instead of hard threshold
% 
% nTrials = length(contrast_ratios);
% 
% out = predict_tRightResponse;
% right_thresh = zeros(1,nTrials);
% labels = zeros(1,nTrials);
% for i = 1:nTrials
%   labels(i) = rand(1);
%   if labels(i)<out(i)
%     right_thresh(i) = 1;
%   end
% end
% 
% right_thresh = right_thresh';
% predict_parameter_table = table(right_thresh, tRatios, tContrasts, tRatios_t_tContrasts);
% 
%  
%  % compare test and predict (% right)
% 
% nActualRightResponseB2 = [];
% totalTrials = [];
% for i = 1:length(unique_contrast_ratios)
%     nActualRightResponseB2(i) = sum(tRightResponse(contrast_ratios == unique_contrast_ratios(i))); %correct right responses per stimulus level
%     totalTrials(i) = sum(contrast_ratios == unique_contrast_ratios(i)); %total right trials per stimulus level
% end
% 
% nPredictRightResponseB2 = [];
% 
% for i = 1:length(unique_contrast_ratios)
%     nPredictRightResponseB2(i) = sum(right_thresh(contrast_ratios == unique_contrast_ratios(i))); % predicted right responses per stimulus level
% end
% 
% 
%  %% plot B2 data      
% 
% %find success per contrast ratio and totals for binomial fit
% 
% 
% [phat, pci] = binofit(nActualRightResponseB2, totalTrials);
% 
% neg = [];
% pos =[];
% 
% for i= 1:length(phat)
%     neg(i) = phat(i) - pci(i,1);
%     pos(i) = pci(i,2) - phat(i);
% end
% 
% 
% % make plot for B2 test data...(actual)
% figure(2)
% % subplot(1,2,1);
% 
% B2TestPlot = errorbar(unique_contrast_ratios, phat, neg, pos, 'Linewidth', 1.5, 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', 'Test');
% hold on 
% % set(gca,'xscale','log')
% % vH = plot([0,0],[0,1]);
% % set(vH, 'Color', 'g');
% % set(gca,'xscale','log')
% % xlabel('Contrast Ratio (R/L)')
% % ylabel('% Right')
% % grid on
% % title('LED % Right by Contrast')
% 
% 
% %%
%  
% %...and with predicted data (B2 model)
%  
% [mphat, mpci] = binofit(nPredictRightResponseB2, totalTrials);
% 
% neg = [];
% pos =[];
% 
% for i= 1:length(mphat)
%     neg(i) = mphat(i) - mpci(i,1);
%     pos(i) = mpci(i,2) - mphat(i);
% end
% 
% figure(2)
% % subplot(1,2,2);
% 
% B2modelPlot = errorbar(unique_contrast_ratios, mphat, neg, pos, 'Linewidth', 1.5, 'Marker', 's', 'MarkerSize', 8, 'Color', 'k', 'DisplayName', 'Testing' );
% hold on
% vH = plot([0,0],[0,1]);
% set(vH, 'Color', 'g');
% % set(gca,'xscale','log')
% xlabel('Contrast Ratio (R/L)')
% ylabel('% Right')
% grid on
% title([mouse '- LED: n = '  num2str(B2nTrials)])
% 
% hold off
% 
% 
% legend([B2TestPlot B2modelPlot], {'Actual', 'Model'}, 'Location', 'northwest')
% 
% % print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
% % mouse '_B2Model_summary.pdf'],'-dpdf','-bestfit');
% % 
% % saveas(gcf,[mouse '_B2Model_summary', '.png'])
% 
% 
% %% Explained Variance
% 
% YB1 = nActualRightResponseB1;
% UB1 = nPredictRightResponseB1;
% 
% YB2 = nActualRightResponseB2;
% UB2 = nPredictRightResponseB2;
% 
% EVB1 = 1-(var(YB1-UB1)/var(YB1))
% EVB2 = 1-(var(YB2-UB2)/var(YB2))
% 
% 
% 
% 
% %% Save workspace
% % 
% % filename = [mouse 'model_workspace.mat'];
% % save(filename)