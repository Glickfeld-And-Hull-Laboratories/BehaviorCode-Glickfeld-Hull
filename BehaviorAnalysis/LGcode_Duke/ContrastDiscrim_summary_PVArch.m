close all
clear all

mouse_list = {'i441','i442','i443'};
nmouse = size(mouse_list,2);
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC';
b2Ix_all = [];
tRight_all = [];
ratio_mat_all = [];
IiX_all = [];
rCon_all = [];
thresh_all = zeros(nmouse,2);
slope_all = zeros(nmouse,2);
for imouse = 1:nmouse
    mouse = mouse_list{imouse};
    load(fullfile(data_path, [mouse '_trialInfo.mat']));
    rCon_all = [rCon_all rCon];
    b2Ix_all = [b2Ix_all b2Ix];
    tRight_all = [tRight_all tRight];
    ratio_mat_all = [ratio_mat_all ratio_mat];
    IiX_all = [IiX_all IiX];
    thresh_all(imouse,:) = thresh_mat;
    slope_all(imouse,:) = slope_mat;
end
thresh_all(find(thresh_all<1)) = -1./thresh_all(find(thresh_all<1));

ratios = unique(ratio_mat_all);
tot_all = zeros(2, length(ratios));
totR_all = zeros(2, length(ratios));
pct_right_all = zeros(2,length(ratios));
ci_right_all = zeros(2,2,length(ratios));
for irat = 1:length(ratios)
    for ib = 1:2
        tot_all(ib,irat) = sum(ratio_mat_all==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1);
        totR_all(ib,irat) = sum(ratio_mat_all==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & tRight_all);
        [pct_right_all(ib,irat) ci_right_all(:,ib,irat)] = binofit(totR_all(ib,irat),tot_all(ib,irat));
    end
end

figure
col_mat = [[0 0 0]; defaultPlotColors(1)];
for ib = 1:2
    scatter(ratios, pct_right_all(ib,:), tot_all(ib,:)./5, col_mat(ib,:))
    hold on
end
set(gca, 'XScale', 'log')
xlabel('Contrast Ratio')
ylabel('Fraction Right Choice')
title(cell2mat(mouse_list))
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' cell2mat(mouse_list) '_Summary.pdf'], '-dpdf','-bestfit')
%%
figure;
%thresh
subplot(1,2,1)
scatter(ones(nmouse,1),thresh_all(:,2)-thresh_all(:,1), 'o')
hold on
errorbar(1,mean(thresh_all(:,2)-thresh_all(:,1),1), std(thresh_all(:,2)-thresh_all(:,1),1)./sqrt(nmouse),'o')
xlim([0 2])
ylim([-4 4])
ylabel('Change in threshold')

%thresh
subplot(1,2,2)
scatter(ones(nmouse,1),slope_all(:,2)-slope_all(:,1), 'o')
hold on
errorbar(1,mean(slope_all(:,2)-slope_all(:,1),1), std(slope_all(:,2)-slope_all(:,1),1)./sqrt(nmouse),'o')
xlim([0 2])
ylim([-2 2])
ylabel('Change in slope')

suptitle(cell2mat(mouse_list))
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' cell2mat(mouse_list) '_Thresh&SlopeSummary.pdf'], '-dpdf','-bestfit')

%%

edges = [0.0095 0.04 0.2 0.5 1 1.25 2.5 5 80 102];
[bin n idx] = histcounts(ratio_mat_all, edges);
tRatBin = zeros(size(ratio_mat_all));
for i = 1:length(edges)
    tRatBin(find(idx==i)) = mean(ratio_mat_all(find(idx==i)));
end

ratios = unique(tRatBin);
tot_all = zeros(2, length(ratios));
totR_all = zeros(2, length(ratios));
pct_right_all = zeros(2,length(ratios));
ci_right_all = zeros(2,2,length(ratios));
for irat = 1:length(ratios)
    for ib = 1:2
        tot_all(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1);
        totR_all(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & tRight_all);
        [pct_right_all(ib,irat) ci_right_all(:,ib,irat)] = binofit(totR_all(ib,irat),tot_all(ib,irat));
    end
end

figure;
col_mat = ['k'; 'c'];
for ib = 1:2
    errorbar(ratios, pct_right_all(ib,:), pct_right_all(ib,:)-squeeze(ci_right_all(1,ib,:))', squeeze(ci_right_all(2,ib,:))'-pct_right_all(ib,:), ['o' col_mat(ib)]) 
    hold on
end
set(gca, 'XScale', 'log')

%% by right con
thresh = 0.4;
ratio_mat_all(find(ratio_mat_all<0.01)) = 0.01;
ratio_mat_all(find(ratio_mat_all>100)) = 100;
edges = [0.0095 0.04 0.2 0.5 1 1.25 2.5 5 80 102];
[bin n idx] = histcounts(ratio_mat_all, edges);
tRatBin = zeros(size(ratio_mat_all));
for i = 1:length(edges)
    tRatBin(find(idx==i)) = mean(ratio_mat_all(find(idx==i)));
end
ratios = unique(tRatBin);
tot_all_high = zeros(2, length(ratios));
totR_all_high = zeros(2, length(ratios));
pct_right_all_high = zeros(2,length(ratios));
ci_right_all_high = zeros(2,2,length(ratios));
tot_all_low = zeros(2, length(ratios));
totR_all_low = zeros(2, length(ratios));
pct_right_all_low = zeros(2,length(ratios));
ci_right_all_low = zeros(2,2,length(ratios));
for irat = 1:length(ratios)
    for ib = 1:2
        tot_all_low(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & rCon_all<thresh);
        tot_all_high(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & rCon_all>=thresh);
        totR_all_low(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & tRight_all & rCon_all<thresh);
        totR_all_high(ib,irat) = sum(tRatBin==ratios(irat) & IiX_all==0 & b2Ix_all==ib-1 & tRight_all & rCon_all>=thresh);
        [pct_right_all_high(ib,irat) ci_right_all_high(:,ib,irat)] = binofit(totR_all_high(ib,irat),tot_all_high(ib,irat));
        [pct_right_all_low(ib,irat) ci_right_all_low(:,ib,irat)] = binofit(totR_all_low(ib,irat),tot_all_low(ib,irat));
    end
end

figure;
col_mat = ['k'; 'c'];
for ib = 1:2
    subplot(1,2,1)
    errorbar(ratios, pct_right_all_high(ib,:), pct_right_all_high(ib,:)-squeeze(ci_right_all_high(1,ib,:))', squeeze(ci_right_all_high(2,ib,:))'-pct_right_all_high(ib,:), ['o' col_mat(ib)]) 
    hold on
    subplot(1,2,2)
    errorbar(ratios, pct_right_all_low(ib,:), pct_right_all_low(ib,:)-squeeze(ci_right_all_low(1,ib,:))', squeeze(ci_right_all_low(2,ib,:))'-pct_right_all_low(ib,:), ['o' col_mat(ib)]) 
    hold on
end
subplot(1,2,1)
set(gca, 'XScale', 'log')
xlabel('Contrast ratio')
ylabel('Fraction right choice')
title('High Contrast')
xlim([0.01 100])
ylim([0 1])
subplot(1,2,2)
set(gca, 'XScale', 'log')
xlabel('Contrast ratio')
ylabel('Fraction right choice')
title('Low Contrast')
xlim([0.01 100])
ylim([0 1])
suptitle(cell2mat(mouse_list))
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' cell2mat(mouse_list) '_HighVsLowCon.pdf'], '-dpdf','-bestfit')
