datapath = 'Z:\home\andrew\Behavior\Data';
cd(datapath)
mouse_list = [508; 509; 511; 512];
s = struct;
for imouse = 1:length(mouse_list)
list = dir(['data-i' num2str(mouse_list(imouse,:)) '-*']);
siz = size(list);
ndays = 0;
pct_corr = [];
pct_miss = [];
pct_early = [];
flash = [];
reqhold = [];
min_ori = [];
for ifile = 1:siz(1)
    load(fullfile(datapath,list(ifile).name));
    if isfield(input, 'reactTimeMs')
        ndays = ndays+1;
        A=cell2mat(input.reactTimesMs);
        ind_corr = intersect(find(A>100),find(A<550));
        ind_miss = find(A>550);
        ind_early = find(A<100);
        pct_corr = [pct_corr length(ind_corr)./length(A)];
        pct_miss = [pct_miss length(ind_miss)./length(A)];
        pct_early = [pct_early length(ind_early)./length(A)];
        if isfield(input, 'stimOnTimeMs')
            flash = [flash 1];
        else
            flash = [flash 0];
        end
        if isfield(input, 'randReqHoldMaxMs')
            reqhold = [reqhold (input.randReqHoldMaxMs + input.fixedReqHoldTimeMs)];
        else
            reqhold = [reqhold ((input.stimOnTimeMs+input.stimOffTimeMs)*input.maxCyclesOn)];
        end
        min_ori = [min_ori min(cell2mat(input.tGratingDirectionDeg),[],2)];
    end
end

figure;
subplot(3,1,1)
plot(pct_corr, 'k');
hold on
plot(pct_miss, 'r');
hold on
plot(pct_early, 'c');
hold on
plot(1+(flash.*0.1), 'b')
subplot(3,1,2)
plot(reqhold);
subplot(3,1,3)
plot(min_ori);
suptitle(num2str(mouse_list(imouse,:)))

s(imouse).pct_corr = pct_corr;
s(imouse).pct_miss = pct_miss;
s(imouse).pct_early = pct_early;
s(imouse).flash = flash;
s(imouse).reqhold = reqhold;
s(imouse).min_ori = min_ori;
end

preflash =[];
postflash = [];
for imouse = 1:length(mouse_list)
    ind_pre = find(s(imouse).flash == 0);
    ind_post = find(s(imouse).flash == 1);
    preflash = [preflash length(ind_pre)];
    postflash = [postflash length(ind_post)];
end
max_pre = max(preflash,[],2);
max_post = max(postflash,[],2);

pct_corr_pre_all = zeros(length(mouse_list), max_pre);
pct_miss_pre_all = zeros(length(mouse_list), max_pre);
pct_early_pre_all = zeros(length(mouse_list), max_pre);
reqhold_pre_all = zeros(length(mouse_list), max_pre);
min_ori_pre_all = zeros(length(mouse_list), max_pre);

pct_corr_post_all = zeros(length(mouse_list), max_post);
pct_miss_post_all = zeros(length(mouse_list), max_post);
pct_early_post_all = zeros(length(mouse_list), max_post);
reqhold_post_all = zeros(length(mouse_list), max_post);
min_ori_post_all = zeros(length(mouse_list), max_post);

for imouse = 1:length(mouse_list)
    ind_pre = find(s(imouse).flash == 0);
    ind_post = find(s(imouse).flash == 1);
    pct_corr_pre_all(imouse,1:length(ind_pre)) = s(imouse).pct_corr(1:length(ind_pre));
    pct_miss_pre_all(imouse,1:length(ind_pre)) = s(imouse).pct_miss(1:length(ind_pre));
    pct_early_pre_all(imouse,1:length(ind_pre)) = s(imouse).pct_early(1:length(ind_pre));
    reqhold_pre_all(imouse,1:length(ind_pre)) = s(imouse).reqhold(1:length(ind_pre));
    min_ori_pre_all(imouse,1:length(ind_pre)) = s(imouse).min_ori(1:length(ind_pre));
    
    pct_corr_post_all(imouse,1:length(ind_post)) = s(imouse).pct_corr(ind_post(1):ind_post(1)+length(ind_post)-1);
    pct_miss_post_all(imouse,1:length(ind_post)) = s(imouse).pct_miss(ind_post(1):ind_post(1)+length(ind_post)-1);
    pct_early_post_all(imouse,1:length(ind_post)) = s(imouse).pct_early(ind_post(1):ind_post(1)+length(ind_post)-1);
    reqhold_post_all(imouse,1:length(ind_post)) = s(imouse).reqhold(ind_post(1):ind_post(1)+length(ind_post)-1);
    min_ori_post_all(imouse,1:length(ind_post)) = s(imouse).min_ori(ind_post(1):ind_post(1)+length(ind_post)-1);
    
    if length(ind_pre<max_pre)
        pct_corr_pre_all(imouse,length(ind_pre)+1:max_pre) = NaN;
        pct_miss_pre_all(imouse,length(ind_pre)+1:max_pre) = NaN;
        pct_early_pre_all(imouse,length(ind_pre)+1:max_pre) = NaN;
        reqhold_pre_all(imouse,length(ind_pre)+1:max_pre) = NaN;
        min_ori_pre_all(imouse,length(ind_pre)+1:max_pre) = NaN;
    end
    if length(ind_post<max_post)
        pct_corr_post_all(imouse,length(ind_post)+1:max_post) = NaN;
        pct_miss_post_all(imouse,length(ind_post)+1:max_post) = NaN;
        pct_early_post_all(imouse,length(ind_post)+1:max_post) = NaN;
        reqhold_post_all(imouse,length(ind_post)+1:max_post) = NaN;
        min_ori_post_all(imouse,length(ind_post)+1:max_post) = NaN;
    end
end
pct_corr_pre_avg = nanmean(pct_corr_pre_all,1);
pct_corr_pre_std = nanstd(pct_corr_pre_all,[],1);
pct_miss_pre_avg = nanmean(pct_miss_pre_all,1);
pct_miss_pre_std = nanstd(pct_miss_pre_all,[],1);
pct_early_pre_avg = nanmean(pct_early_pre_all,1);
pct_early_pre_std = nanstd(pct_early_pre_all,[],1);
reqhold_pre_avg = nanmean(reqhold_pre_all,1);
reqhold_pre_std = nanstd(reqhold_pre_all,[],1);
min_ori_pre_avg = nanmean(min_ori_pre_all,1);
min_ori_pre_std = nanstd(min_ori_pre_all,[],1);

figure;
subplot(2,1,1)
errorbar(pct_corr_pre_avg, pct_corr_pre_std,'k');
hold on
errorbar(pct_miss_pre_avg, pct_miss_pre_std,'r');
hold on
errorbar(pct_early_pre_avg, pct_early_pre_std,'c');
subplot(2,1,2)
errorbar(reqhold_pre_avg, reqhold_pre_std);
suptitle('Early phase')

pct_corr_post_avg = nanmean(pct_corr_post_all,1);
pct_corr_post_std = nanstd(pct_corr_post_all,[],1);
pct_miss_post_avg = nanmean(pct_miss_post_all,1);
pct_miss_post_std = nanstd(pct_miss_post_all,[],1);
pct_early_post_avg = nanmean(pct_early_post_all,1);
pct_early_post_std = nanstd(pct_early_post_all,[],1);
reqhold_post_avg = nanmean(reqhold_post_all,1);
reqhold_post_std = nanstd(reqhold_post_all,[],1);
min_ori_post_avg = nanmean(min_ori_post_all,1);
min_ori_post_std = nanstd(min_ori_post_all,[],1);

figure;
subplot(2,1,1)
errorbar(pct_corr_post_avg, pct_corr_post_std,'k');
hold on
errorbar(pct_miss_post_avg, pct_miss_post_std,'r');
hold on
errorbar(pct_early_post_avg, pct_early_post_std,'c');
subplot(2,1,2)
errorbar(min_ori_post_avg, min_ori_post_std);
suptitle('Late phase')

figure;
col_mat = strvcat('k','r','b','g');
for imouse = 1:length(mouse_list)
    ind_pre = find(s(imouse).flash == 0);
    subplot(4,1,1)
    plot(s(imouse).pct_corr(1,1:length(ind_pre)),col_mat(imouse,:));
    hold on
    subplot(4,1,2)
    plot(s(imouse).pct_miss(1,1:length(ind_pre)),col_mat(imouse,:));
    hold on
    subplot(4,1,3)
    plot(s(imouse).pct_early(1,1:length(ind_pre)),col_mat(imouse,:));
    hold on
    subplot(4,1,4)
    plot(s(imouse).reqhold(1,1:length(ind_pre)),col_mat(imouse,:));
    hold on
end
suptitle('Early phase')

figure;
col_mat = strvcat('k','r','b','g');
for imouse = 1:length(mouse_list)
    ind_post = find(s(imouse).flash == 1);
    subplot(4,1,1)
    plot(s(imouse).pct_corr(1,ind_post(1):ind_post(1)+length(ind_post)-1),col_mat(imouse,:));
    hold on
    subplot(4,1,2)
    plot(s(imouse).pct_miss(1,ind_post(1):ind_post(1)+length(ind_post)-1),col_mat(imouse,:));
    hold on
    subplot(4,1,3)
    plot(s(imouse).pct_early(1,ind_post(1):ind_post(1)+length(ind_post)-1),col_mat(imouse,:));
    hold on
    subplot(4,1,4)
    plot(s(imouse).min_ori(1,ind_post(1):ind_post(1)+length(ind_post)-1),col_mat(imouse,:));
    hold on
end
suptitle('Late phase')