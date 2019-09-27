%Begin with input_for_pSignifit, get log fit & slopes 
%Plot individually then find avg of all animals (include error bars)

%Get decision thresholds (50% correct or left)
%result, result0, and result90 are psignifit outputs 
%pCorrect is 0.5
%don't understand unscaled, try both
pCorrect= 0.5;
[thresh_NA_1401,CI_NA_1401] = getThreshold(result, pCorrect, 'unscaled')
[thresh_90_1401,CI_90_1401] = getThreshold(result90, pCorrect, 'unscaled')
[thresh_0_1401,CI_0_1401] = getThreshold(result0, pCorrect, 'unscaled')


% [thresh_NA_1402,CI_NA_1402] = getThreshold(result, pCorrect, 'unscaled')
% [thresh_90_1402,CI_90_1402] = getThreshold(result90, pCorrect, 'unscaled')
% [thresh_0_1402,CI_0_1402] = getThreshold(result0, pCorrect, 'unscaled')


%Save DecisionThresholds & psych curve slope @ dec thresh for each animal.
slope_max_NA_1401 = slope_max;
slope_max_90_1401 = slope_max90;
slope_max_0_1401 = slope_max0;

slope_max_NA_1402 = slope_max;
slope_max_90_1402 = slope_max90;
slope_max_0_1402 = slope_max0;


%%Quantify difference in max slope
conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});

% max_slopes_1401 = horzcat(slope_max_NA_1401,slope_max_90_1401,slope_max_0_1401);
dec_thresh_1401 = horzcat(thresh_NA_1401, thresh_90_1401, thresh_0_1401);
CI_thresh_1401 = vertcat(CI_NA_1402, CI_90_1402, CI_0_1402)';


% max_slopes_1402 = horzcat(slope_max_NA_1402,slope_max_90_1402,slope_max_0_1402);
% dec_thresh_1402 = horzcat(thresh_NA_1402, thresh_90_1402, thresh_0_1402);
% CI_thresh_1402 = vertcat(CI_NA_1402, CI_90_1402, CI_0_1402)';


% avg_max_slopes = (max_slopes_1401 + max_slopes_1402)./2
% max_slopes_NA = horzcat(max_slopes_1401(1),max_slopes_1402(1));
% max_slopes_90 = horzcat(max_slopes_1401(2),max_slopes_1402(2));
% max_slopes_0 = horzcat(max_slopes_1401(3),max_slopes_1402(3));
% std_max_slopes_NA = std(max_slopes_NA)/sqrt(2);
% std_max_slopes_90 = std(max_slopes_90)/sqrt(2);
% std_max_slopes_0 = std(max_slopes_0)/sqrt(2);
% stds_max_slopes = horzcat(std_max_slopes_NA, std_max_slopes_90,std_max_slopes_0);

% avg_dec_thresh = (dec_thresh_1401 + dec_thresh_1402)./2
% dec_thresh_NA = horzcat(dec_thresh_1401(1),dec_thresh_1402(1));
% dec_thresh_90 = horzcat(dec_thresh_1401(2),dec_thresh_1402(2));
% dec_thresh_0 = horzcat(dec_thresh_1401(3),dec_thresh_1402(3));
% std_dec_thresh_NA = std(dec_thresh_NA)/sqrt(2);
% std_dec_thresh_90 = std(dec_thresh_90)/sqrt(2);
% std_dec_thresh_0 = std(dec_thresh_0)/sqrt(2);
% stds_dec_thresh = horzcat(std_dec_thresh_NA,std_dec_thresh_90,std_dec_thresh_0);



% figure
% plot(conditions,max_slopes_1401,'-o')
% hold on
% plot(conditions,max_slopes_1402,'-o')
% errorbar(conditions,avg_max_slopes,stds_max_slopes)
% legend('i1401','i1402','Avg +/- SEM')
% title({'Max Slopes of Psychometric Fits Across Conditions'; 'Subjects: i1401,i1402'})
% xlabel('Adaptation Condition')
% ylabel('Max Slope of Psych Curve Fit')
% hold off

% figure
% plot(conditions,dec_thresh_1401,'-o')
% hold on
% plot(conditions,dec_thresh_1402,'-o')
% legend('i1401','i1402')
% title({'Discrimination Threshold based on Psychometric Fits Across Conditions';'Subjects: i1401,i1402'})
% xlabel('Adaptation Condition')
% ylabel('Discrimination Threshold')
% ylim([5,13])
% hold off

figure
plot(conditions, dec_thresh_1401,'-o','Color',[0 0.4470 0.7410])
hold on
f = errorbar(conditions,dec_thresh_1401,(dec_thresh_1402-CI_thresh_1401(1,:)),(CI_thresh_1401(2,:)-dec_thresh_1402),'-o')
f.Color = [0 0.4470 0.7410];
plot(conditions,dec_thresh_1402,'-o','Color',[0.8500 0.3250 0.0980])
e= errorbar(conditions,dec_thresh_1402,(dec_thresh_1402-CI_thresh_1402(1,:)),(CI_thresh_1402(2,:)-dec_thresh_1402),'-o')
e.Color = [0.8500 0.3250 0.0980];
legend('i1401','i1401','i1402','i1402')
title({'Discrimination Threshold based on Psychometric Fits Across Conditions';'Subjects: i1401,i1402'})
xlabel('Adaptation Condition')
ylabel('Discrimination Threshold')
ylim([0,20])
hold off

%% Alternate approach to calculating threshold - ori at peak slope 

% alt_thresh_NA_1402 = Stimlevel(max_idx);
% alt_thresh_90_1402 = Stimlevel(max_idx90);
% alt_thresh_0_1402 = Stimlevel(max_idx0);
% alt_thresh_1402 = horzcat(alt_thresh_NA_1402, alt_thresh_90_1402, alt_thresh_0_1402);

alt_thresh_NA_1401 = Stimlevel(max_idx);
alt_thresh_90_1401 = Stimlevel(max_idx90);
alt_thresh_0_1401 = Stimlevel(max_idx0);
alt_thresh_1401 = horzcat(alt_thresh_NA_1401, alt_thresh_90_1401, alt_thresh_0_1401);


avg_alt_thresh = (alt_thresh_1401 + alt_thresh_1402)./2
alt_thresh_NA = horzcat(alt_thresh_1401(1),alt_thresh_1402(1));
alt_thresh_90 = horzcat(alt_thresh_1401(2),alt_thresh_1402(2));
alt_thresh_0 = horzcat(alt_thresh_1401(3),alt_thresh_1402(3));
std_alt_thresh_NA = std(alt_thresh_NA)/sqrt(2);
std_alt_thresh_90 = std(alt_thresh_90)/sqrt(2);
std_alt_thresh_0 = std(alt_thresh_0)/sqrt(2);
stds_alt_thresh = horzcat(std_alt_thresh_NA,std_alt_thresh_90,std_alt_thresh_0);

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});

figure
plot(conditions,alt_thresh_1401,'-o')
hold on
plot(conditions,alt_thresh_1402,'-o')
errorbar(conditions,avg_alt_thresh,stds_alt_thresh)
legend('i1401','i1402','Avg +/- SEM')
title({'Alternative Decision Threshold based on Psychometric Fits Across Conditions';'Subjects: i1401,i1402'})
xlabel('Adaptation Condition')
ylabel('Decision Threshold')
hold off

%% Delta change in threshold in 90 & 0 deg adapter conditions as compared to NA condition - use standard thrseh (50% correct) 

% dec_thresh_1402 = horzcat(thresh_NA_1402, thresh_90_1402, thresh_0_1402);
% dec_thresh_1401 = horzcat(thresh_NA_1401, thresh_90_1401, thresh_0_1401);

delta_thresh_1401 = [abs(dec_thresh_1401(2))/abs(dec_thresh_1401(1)) abs(dec_thresh_1401(3))/abs(dec_thresh_1401(1))];
delta_thresh_1402 = [abs(dec_thresh_1402(2))/abs(dec_thresh_1402(1)) abs(dec_thresh_1402(3))/abs(dec_thresh_1402(1))];

conditions_90_0 = categorical(["90deg Adapter","0deg Adapter"]);
conditions_90_0 = reordercats(conditions_90_0,{'90deg Adapter','0deg Adapter'});

figure
plot(conditions_90_0,delta_thresh_1401,'-o')
hold on
plot(conditions_90_0,delta_thresh_1402,'-o')
legend('i1401','i1402')
title({'Absolute Value of Delta Change in Decision Threshold in Presence of Adapter'})
xlabel('Adaptation Condition')
ylabel('Absolute Delta Change Decision Threshold')
hold off 

%% Delta change in max slope in 90 & 0 deg adapter conditions
% I think CI's are wong

%max_slopes_1401 = horzcat(slope_max_NA_1401,slope_max_90_1401,slope_max_0_1401);
%max_slopes_1402 = horzcat(slope_max_NA_1402,slope_max_90_1402,slope_max_0_1402);

conditions_90_0 = categorical(["90deg Adapter","0deg Adapter"]);
conditions_90_0 = reordercats(conditions_90_0,{'90deg Adapter','0deg Adapter'});

delta_slope_1401 = [max_slopes_1401(2)/max_slopes_1401(1) max_slopes_1401(3)/max_slopes_1401(1)]
delta_slope_1402 = [max_slopes_1402(2)/max_slopes_1402(1) max_slopes_1402(3)/max_slopes_1402(1)];

%1401 CIs for delta slope (using bootstrap) - I think this is wrong; have
%to redo whole bootstrap to get CIs.... 
slopies_1401_NA = sort(slopies_1401(:,1));
slopies_1401_NA_5 = slopies_1401_NA(5)
slopies_1401_NA_95 = slopies_1401_NA(95)

slopies_1401_90 = sort(slopies_1401(:,2));
slopies_1401_90_5 = slopies_1401_90(5)
slopies_1401_90_95 = slopies_1401_90(95)
delta_slopies_1401_90_5 = slopies_1401_90_5/slopies_1401_NA_5
delta_slopies_1401_90_95 = slopies_1401_90_95/slopies_1401_NA_5
CI_1401_90 = vertcat(delta_slopies_1401_90_5,delta_slopies_1401_90_95)

slopies_1401_0 = sort(slopies_1401(:,3));
slopies_1401_0_5 = slopies_1401_0(5)
slopies_1401_0_95 = slopies_1401_0(95)
delta_slopies_1401_0_5 = slopies_1401_0_5/slopies_1401_NA_5
delta_slopies_1401_0_95 = slopies_1401_0_95/slopies_1401_NA_5
CI_1401_0 = vertcat(delta_slopies_1401_0_5,delta_slopies_1401_0_95)

CI_1401_90_0 = horzcat(CI_1401_90,CI_1401_0);

%1402 CIs for delta slopes
slopies_1402_NA = sort(slopies_1402(:,1));
slopies_1402_NA_5 = slopies_1402_NA(5)
slopies_1402_NA_95 = slopies_1402_NA(95)

slopies_1402_90 = sort(slopies_1402(:,2));
slopies_1402_90_5 = slopies_1402_90(5)
slopies_1402_90_95 = slopies_1402_90(95)
delta_slopies_1402_90_5 = slopies_1402_90_5/slopies_1402_NA_5
delta_slopies_1402_90_95 = slopies_1402_90_95/slopies_1402_NA_5
CI_1402_90 = vertcat(delta_slopies_1402_90_5,delta_slopies_1402_90_95)

slopies_1402_0 = sort(slopies_1402(:,3));
slopies_1402_0_5 = slopies_1402_0(5)
slopies_1402_0_95 = slopies_1402_0(95)
delta_slopies_1402_0_5 = slopies_1402_0_5/slopies_1402_NA_5
delta_slopies_1402_0_95 = slopies_1402_0_95/slopies_1402_NA_5
CI_1402_0 = vertcat(delta_slopies_1402_0_5,delta_slopies_1402_0_95)

CI_1402_90_0 = horzcat(CI_1402_90, CI_1402_0);

figure
plot(conditions_90_0,delta_slope_1401,'-o','Color',[0 0.4470 0.7410])
hold on
% f = errorbar(conditions_90_0,delta_slope_1401,(delta_slope_1401-CI_1401_90_0(1,:)),(CI_1401_90_0(2,:)-delta_slope_1401),'-o')
f = errorbar(conditions_90_0,delta_slope_1401,(CI_1401_90_0(:,1)-delta_slope_1401(1)),(CI_1401_90_0(:,2)-delta_slope_1401(2)),'-o')
f.Color = [0 0.4470 0.7410];
plot(conditions_90_0,delta_slope_1402,'-o')
e= errorbar(conditions_90_0, delta_slope_1402,(CI_1402_90_0(:,1)-delta_slope_1402(1)),(CI_1402_90_0(:,2)-delta_slope_1402(2)),'-o')
e.Color = [0.8500 0.3250 0.0980];
legend('i1401','i1401','i1402','i1402')
title({'Delta Change in Max Slope in Presence of Adapter'})
xlabel('Adaptation Condition')
ylabel('Delta Change Decision Threshold')
ylim([-0.5,2])
hold off 






