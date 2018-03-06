mouse_str = {'519','520','523','524','525','528','529','530','531','532','535','537','538','628','547','548','550','551','552','553','560','561','564','565','566','567','575','578','579','580','582','583','586'};
nmouse = size(mouse_str,2);

for imouse = 1:nmouse
    fn = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\';
    mouse = mouse_str{imouse};
    allDays = dir([fn 'data-i' mouse '-*']);
    ndates = size(allDays,1);
    dates = [];
    for iday = 1:ndates
        dates = [dates; str2num(allDays(iday).name(1,11:16))];
    end
    udates = unique(dates);
    nudates = length(udates);
    dayIx = struct;
    ndays = 0;
    for i = 1:nudates
        if ndays<31
            fprintf([num2str(ndays) ' '])
            ind = find(dates == udates(i));
            if length(ind)>1
                ds = struct;
                for ii = 1:length(ind)
                    load([fn allDays(ind(ii)).name]);
                    if isfield(input, 'tLeftTrial')
                        if size(input.trialOutcomeCell,1)>20
                            ds(ii) = input; 
                        end
                    end
                end
                if isfield(ds, 'tLeftTrial')
                    input = concatenateDataBlocks(ds);
                end
            else
                load([fn allDays(ind).name]);
            end
            if isfield(input, 'tLeftTrial')
                correctIx = strcmp(input.trialOutcomeCell, 'success');
                incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
                if isfield(input,'doNoGo')
                    if input.doNoGo
                        isNoGoIx = celleqel2mat_padded(input.isNoGo);
                    else
                        isNoGoIx = zeros(size(correctIx));
                    end
                else
                    isNoGoIx = zeros(size(correctIx));
                end
                corrIx = setdiff(find(correctIx),find(isNoGoIx));
                incorrIx = setdiff(find(incorrectIx),find(isNoGoIx));
                respIx = setdiff(find(correctIx+incorrectIx),find(isNoGoIx));
                if length(respIx)>20
                    ndays = ndays+1;
                    dayIx(ndays).date = udates(i);
                    dayIx(ndays).trialsCompleted = sum(correctIx) + sum(incorrectIx);
                    decisionTimeMs = celleqel2mat_padded(input.tDecisionTimeMs);
                    earlyDecisionTimeMs = mean(decisionTimeMs(1,respIx(1,1:20)),2);
                    if length(respIx)>24
                        for iresp = 21:length(respIx)-4
                            if size(respIx,2)>iresp+4
                                ind = find(decisionTimeMs(1,respIx(1,iresp:iresp+4))>2*earlyDecisionTimeMs);
                                if sum(ind) == 5
                                    endTrial = iresp;
                                end
                            else 
                                endTrial = length(respIx);
                            end
                        end
                    else
                        endTrial = length(respIx);
                    end
                    dayIx(ndays).avgDecisionTimeMs = mean(decisionTimeMs(1,respIx(1,1:endTrial)),2);
                    tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
                    tLeftResponse = zeros(size(tLeftTrial));
                    tLeftResponse(intersect(find(tLeftTrial),find(corrIx))) = 1;
                    tLeftResponse(intersect(find(~tLeftTrial),find(incorrIx))) = 1;
                    tGratingContrast = cell2mat_padded(input.tGratingContrast);
                    nRightRespRightTrial = length(intersect(find(tGratingContrast(1,respIx(1,1:endTrial)) == 1),intersect(find(~tLeftResponse(1,respIx(1,1:endTrial))),find(~tLeftTrial(1,respIx(1,1:endTrial))))));
                    nRightTrial = length(intersect(find(tGratingContrast(1,respIx(1,1:endTrial)) == 1),find(~tLeftTrial(1,respIx(1,1:endTrial)))));
                    nRightRespLeftTrial = length(intersect(find(tGratingContrast(1,respIx(1,1:endTrial)) == 1),intersect(find(~tLeftResponse(1,respIx(1,1:endTrial))),find(tLeftTrial(1,respIx(1,1:endTrial))))));
                    nLeftTrial = length(intersect(find(tGratingContrast(1,respIx(1,1:endTrial)) == 1),find(tLeftTrial(1,respIx(1,1:endTrial)))));
                    dayIx(ndays).selectivity = (nRightRespRightTrial/nRightTrial)-(nRightRespLeftTrial/nLeftTrial);
                    dayIx(ndays).bias = abs((((nRightRespRightTrial/nRightTrial)+(nRightRespLeftTrial/nLeftTrial))/2)-0.5);
                end
            end
        end
    end
    fprintf(' \n')
    fn_out = ['\\crash.dhe.duke.edu\data\home\lindsey\Analysis\Behavior\2AFC\i' mouse '_dayIx.mat'];
    save(fn_out,'dayIx');
end

%%
figure;
mouseIx = struct;
blue_mat = blues(100);
trialsCompleted = nan(nmouse,31);
selectivity = nan(nmouse,31);
bias = nan(nmouse,31);
avgDecisionTimeMs = nan(nmouse,31);
nDays = zeros(1,nmouse);
for imouse = 1:nmouse
    mouse = mouse_str{imouse};
    fn_out = ['\\crash.dhe.duke.edu\data\home\lindsey\Analysis\Behavior\2AFC\i' mouse '_dayIx.mat'];
    load(fn_out)
    nDays(1,imouse) = size(dayIx,2);
    subplot(2,2,1)
    trialsCompleted(imouse,1:nDays(1,imouse)) = [dayIx.trialsCompleted];
    plot([dayIx.trialsCompleted],'Color',blue_mat(65+imouse,:))
    hold on
    subplot(2,2,2)
    selectivity(imouse,1:nDays(1,imouse)) = [dayIx.selectivity];
    plot([dayIx.selectivity],'Color',blue_mat(65+imouse,:))
    hold on
    subplot(2,2,3)
    bias(imouse,1:nDays(1,imouse)) = [dayIx.bias];
    plot([dayIx.bias],'Color',blue_mat(65+imouse,:))
    hold on
    subplot(2,2,4)
    avgDecisionTimeMs(imouse,1:nDays(1,imouse)) = [dayIx.avgDecisionTimeMs];
    plot([dayIx.avgDecisionTimeMs],'Color',blue_mat(65+imouse,:))
    hold on
    
end
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 inf])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.5])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 15000])
xlim([0 35])
suptitle(['All 2AFC mice; n = ' num2str(nmouse)])

nanmouse = ~isnan(trialsCompleted);
figure;
subplot(2,2,1)
shadedErrorBar(1:31,nanmean(trialsCompleted,1), nanstd(trialsCompleted,[],1)./(sqrt(sum(nanmouse,1))));
subplot(2,2,2)
shadedErrorBar(1:31,nanmean(selectivity,1), nanstd(selectivity,[],1)./(sqrt(sum(nanmouse,1))));
subplot(2,2,3)
shadedErrorBar(1:31,nanmean(bias,1), nanstd(bias,[],1)./(sqrt(sum(nanmouse,1))));
subplot(2,2,4)
shadedErrorBar(1:31,nanmean(avgDecisionTimeMs,1), nanstd(avgDecisionTimeMs,[],1)./(sqrt(sum(nanmouse,1))));
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 350])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-.1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.2])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 10000])
xlim([0 35])
suptitle(['All 2AFC mice; n = ' num2str(nmouse)])

figure;
ind = find(nDays>30);
for imouse = 1:nmouse
    if find(ind==imouse)
        subplot(2,2,1)
        plot(trialsCompleted(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,2)
        plot(selectivity(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,3)
        plot(bias(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,4)
        plot(avgDecisionTimeMs(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
    end
end
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 inf])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.5])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 15000])
xlim([0 35])
suptitle(['All 2AFC mice with 31 days training; n = ' num2str(length(ind))])

figure;
subplot(2,2,1)
shadedErrorBar(1:31,nanmean(trialsCompleted(ind,:),1), nanstd(trialsCompleted(ind,:),[],1)./(sqrt(sum(nanmouse(ind,:),1))));
subplot(2,2,2)
shadedErrorBar(1:31,nanmean(selectivity(ind,:),1), nanstd(selectivity(ind,:),[],1)./(sqrt(sum(nanmouse(ind,:),1))));
subplot(2,2,3)
shadedErrorBar(1:31,nanmean(bias(ind,:),1), nanstd(bias(ind,:),[],1)./(sqrt(sum(nanmouse(ind,:),1))));
subplot(2,2,4)
shadedErrorBar(1:31,nanmean(avgDecisionTimeMs(ind,:),1), nanstd(avgDecisionTimeMs(ind,:),[],1)./(sqrt(sum(nanmouse(ind,:),1))));
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 350])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-.1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.2])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 10000])
xlim([0 35])
suptitle(['All 2AFC mice with 31 days training; n = ' num2str(length(ind))])

mouse_use = [1,0,0,1,0,0,0,1,1,0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0];
mouse_ind = find(mouse_use);
nmouse_used = sum(mouse_use,2);
figure;
for imouse = 1:nmouse
    if mouse_use(1,imouse)
        subplot(2,2,1)
        plot(trialsCompleted(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,2)
        plot(selectivity(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,3)
        plot(bias(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,4)
        plot(avgDecisionTimeMs(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
    end
end
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 inf])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.5])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 15000])
xlim([0 35])
suptitle(['All 2AFC mice that learned; n = ' num2str(nmouse_used)])

figure;
subplot(2,2,1)
shadedErrorBar(1:31,nanmean(trialsCompleted(mouse_ind,:),1), nanstd(trialsCompleted(mouse_ind,:),[],1)./(sqrt(sum(nanmouse(mouse_ind,:),1))));
subplot(2,2,2)
shadedErrorBar(1:31,nanmean(selectivity(mouse_ind,:),1), nanstd(selectivity(mouse_ind,:),[],1)./(sqrt(sum(nanmouse(mouse_ind,:),1))));
subplot(2,2,3)
shadedErrorBar(1:31,nanmean(bias(mouse_ind,:),1), nanstd(bias(mouse_ind,:),[],1)./(sqrt(sum(nanmouse(mouse_ind,:),1))));
subplot(2,2,4)
shadedErrorBar(1:31,nanmean(avgDecisionTimeMs(mouse_ind,:),1), nanstd(avgDecisionTimeMs(mouse_ind,:),[],1)./(sqrt(sum(nanmouse(mouse_ind,:),1))));
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 350])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-.1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.2])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 10000])
xlim([0 35])
suptitle(['All 2AFC mice that learned; n = ' num2str(nmouse_used)])

mouse_rec = intersect(mouse_ind,17:31);
nmouse_rec = size(mouse_rec,2);
figure;
for imouse = 1:nmouse
    if find(mouse_rec == imouse)
        subplot(2,2,1)
        plot(trialsCompleted(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,2)
        plot(selectivity(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,3)
        plot(bias(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
        subplot(2,2,4)
        plot(avgDecisionTimeMs(imouse,1:nDays(1,imouse)),'Color',blue_mat(65+imouse,:))
        hold on
    end
end
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 inf])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.5])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 15000])
xlim([0 35])
suptitle(['All 2AFC mice that learned; n = ' num2str(nmouse_used)])

figure;
subplot(2,2,1)
shadedErrorBar(1:31,nanmean(trialsCompleted(mouse_rec,:),1), nanstd(trialsCompleted(mouse_rec,:),[],1)./(sqrt(sum(nanmouse(mouse_rec,:),1))));
subplot(2,2,2)
shadedErrorBar(1:31,nanmean(selectivity(mouse_rec,:),1), nanstd(selectivity(mouse_rec,:),[],1)./(sqrt(sum(nanmouse(mouse_rec,:),1))));
subplot(2,2,3)
shadedErrorBar(1:31,nanmean(bias(mouse_rec,:),1), nanstd(bias(mouse_rec,:),[],1)./(sqrt(sum(nanmouse(mouse_rec,:),1))));
subplot(2,2,4)
shadedErrorBar(1:31,nanmean(avgDecisionTimeMs(mouse_rec,:),1), nanstd(avgDecisionTimeMs(mouse_rec,:),[],1)./(sqrt(sum(nanmouse(mouse_rec,:),1))));
subplot(2,2,1)
xlabel('Training day')
ylabel('Number of trials completed')
ylim([0 350])
xlim([0 35])
subplot(2,2,2)
xlabel('Training day')
ylabel('Selectivity')
ylim([-.1 1])
xlim([0 35])
subplot(2,2,3)
xlabel('Training day')
ylabel('Bias')
ylim([0 0.2])
xlim([0 35])
subplot(2,2,4)
xlabel('Training day')
ylabel('Decision time (ms)')
ylim([0 10000])
xlim([0 35])
suptitle(['All 2AFC mice since 8/1/16; n = ' num2str(nmouse_used)])
