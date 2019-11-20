mouse_list = strvcat('i441','i442', 'i443');
nmouse = size(mouse_list,1);
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC';
b2Ix_all = [];
tRight_all = [];
ratio_mat_all = [];
IiX_all = [];
for imouse = 1:nmouse
    mouse = mouse_list(imouse,:);
    load(fullfile(data_path, [mouse '_trialInfo.mat']));
    b2Ix_all = [b2Ix_all b2Ix];
    tRight_all = [tRight_all tRight];
    ratio_mat_all = [ratio_mat_all ratio_mat];
    IiX_all = [IiX_all IiX];
end

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
    scatter(ratios, pct_right_all(ib,:), tot_all(ib,:), col_mat(ib,:))
    hold on
end
set(gca, 'XScale', 'log')

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