clear all
close all

%%
behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
mouse = 'i488';
ds= [mouse '_exptList'];
eval(ds)
doFit = 1;
fprintf([mouse '\n']);

%%
clear temp;
good_exp = [];
for id = 1:size(dates,1)
    fprintf([num2str(dates(id,:))])
    expt_mat = dir(fullfile(behav_path, ['data-' mouse '-' num2str(dates(id,:)) '-*']));
    if size(expt_mat,1) > 1
        if ~isempty(mat_use{id})
            load(fullfile(behav_path,expt_mat(mat_use{id}).name))
        else 
            x = zeros(1,size(expt_mat,1));
            for im = 1:size(expt_mat,1)
                load(fullfile(behav_path,expt_mat(im).name))
                if isfield(input,'trialOutcomeCell')
                    x(1,im) = 1;
                end
            end
            if sum(x,2) == 1
                load(fullfile(behav_path,expt_mat(find(x==1)).name))
            elseif sum(x,2)>1
                error('Too many mat files')
            elseif sum(x,2)==0
                error('No mat files')
            end
        end
    else
        load(fullfile(behav_path,expt_mat.name))
    end
    if ~input.doFeedbackMotion
        if isfield(input,'doEasyStartWithFeedback')
            if input.doEasyStartWithFeedback
                if ~isempty(trials{id})
                    if trials{id}(1)<21
                        trials{id}(1) = 21;
                    end
                else
                    trials{id} = [21 length(input.tGratingContrast)];
                end
            end
        end
        [s b] = selectCalc(input,trials{id});
        fprintf(['- s = ' num2str(chop(s,2)) '; b = ' num2str(chop(b,2))])
        if s>=0.9 & abs(b)<0.1 % s = 0.9 & b < 0.1 is optimal
            if ~isempty(trials{id})
                input = trialChopper(input,trials{id});
            end
            if input.doContrastDiscrim == 1
                temp(id).rightGratingContrast = input.rightGratingContrast;
                temp(id).leftGratingContrast = input.leftGratingContrast;
                
            elseif input.doOriDiscrim == 1
                temp(id).tGratingDirectionStart = input.tGratingDirectionStart;
                
            end
           
            temp(id).tRightResponse = input.tRightResponse;
            temp(id).tBlock2TrialNumber = input.tBlock2TrialNumber;
            temp(id).trialOutcomeCell = input.trialOutcomeCell;
            temp(id).tGratingContrast = input.tGratingContrast;
            temp(id).tTrialLaserPowerMw = double(input.block2TrialLaserPowerMw).*ones(size(celleqel2mat_padded(input.tTrialLaserPowerMw)));
            temp(id).date = dates(id,:).*ones(size(celleqel2mat_padded(input.tTrialLaserPowerMw)));
            temp(id).tLeftTrial = input.tLeftTrial;
            temp(id).tDecisionTimeMs = input.tDecisionTimeMs;
            temp(id).gratingEccentricityDeg = input.gratingEccentricityDeg;
            temp(id).feedbackMotionSensitivity = input.feedbackMotionSensitivity;
            temp(id).tThisTrialStartTimeMs = input.tThisTrialStartTimeMs;
            temp(id).stimTimestampMs = input.stimTimestampMs;
            temp(id).quadratureTimesUs = input.quadratureTimesUs;
            temp(id).block2TrialLaserPowerMw = input.block2TrialLaserPowerMw;
            temp(id).quadratureValues = input.quadratureValues;
            temp(id).delayTimeMs = input.delayTimeMs;
            good_exp = [good_exp id];
            fprintf('- good\n')
        else
            fprintf('\n')
        end
    else
        fprintf('\n')
    end
end


if input.doContrastDiscrim == 1% Contrast
    input = concatenateDataBlocks(temp);    

    rCon = celleqel2mat_padded(input.rightGratingContrast);
    lCon = celleqel2mat_padded(input.leftGratingContrast);
    ratio_mat = chop(rCon./lCon,3);

    % Ratio processing
    if length(unique(ratio_mat)) > 8  
        ratio_mat = wMeanRatios(ratio_mat);
        ratios = unique(ratio_mat);
    end

    ratios = unique(ratio_mat);
   
    tRight = celleqel2mat_padded(input.tRightResponse);
    b2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
    IiX = strcmp(input.trialOutcomeCell,'ignore');
    pow_mat = input.tTrialLaserPowerMw;
    pows = unique(pow_mat);
    nPow = length(pows);
    
    totR = zeros(2, length(ratios),nPow);
    pct_right = zeros(2,length(ratios),nPow);
    ci_right = zeros(2,2,length(ratios),nPow);
    datelist = cell(1,nPow);

    for ipow = 1:nPow
        for irat = 1:length(ratios)
            for ib = 1:2
                totR(ib,irat,ipow) = sum(pow_mat==pows(ipow) & ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1);
                [pct_right(ib,irat,ipow) ci_right(:,ib,irat,ipow)] =...
                    binofit(sum(pow_mat==pows(ipow) & ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight),totR(ib,irat,ipow));
            end
        end
        datelist{ipow} = unique(input.date(find(pow_mat==pows(ipow))));
    end

elseif input.doOriDiscrim == 1 % Ori
    input = concatenateDataBlocks(temp);    

    dirs_mat = celleqel2mat_padded(input.tGratingDirectionStart);
    dirs = unique(dirs_mat);
    

    
    tRight = celleqel2mat_padded(input.tRightResponse);
    b2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
    IiX = strcmp(input.trialOutcomeCell,'ignore');
    pow_mat = input.tTrialLaserPowerMw;
    pows = unique(pow_mat);
    nPow = length(pows);
    
    totR = zeros(2, length(dirs),nPow);
    pct_right = zeros(2,length(dirs),nPow);
    ci_right = zeros(2,2,length(dirs),nPow);
    datelist = cell(1,nPow);
    
    for ipow = 1:nPow
        for idir = 1:length(dirs)
            for ib = 1:2
                totR(ib,idir,ipow) = sum(pow_mat==pows(ipow) & dirs_mat==dirs(idir) & IiX==0 & b2Ix==ib-1);
                [pct_right(ib,idir,ipow) ci_right(:,ib,idir,ipow)] =...
                    binofit(sum(pow_mat==pows(ipow) & dirs_mat==dirs(idir) & IiX==0 & b2Ix==ib-1 & tRight),totR(ib,idir,ipow));
            end
        end
        datelist{ipow} = unique(input.date(find(pow_mat==pows(ipow))));
    end


    
end

n_goodDays = length(good_exp);
sumTrials = length(input.trialOutcomeCell);

% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\'...
%      mouse 'goodData.mat'],'input');

% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Orientation_discrim\'...
%      mouse 'goodData.mat'],'input');

% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\'...
%      mouse '_variables.mat'],'dates', 'trials', 'good_exp', 'mouse', 'pct_right', 'ci_right', 'pows', 'ratios', 'sumTrials', 'n_goodDays',...
%      );


%% plot

% if isfield(input, 'rightGratingContrast') == 1 
%     
%     condition_mat = ratio_mat;
%     u_conditions = ratios;
%     conlabel = 'Contrast ratio (R/L)';
%     
% elseif isfield(input, 'tGratingDirectionStart') == 1
%     
%     condition_mat = dirs_mat;
%     u_conditions = dirs;
%     conlabel = 'Oridiff';
% 
%     
%     
% end

figure;
if nPow > 1
    tPow = nPow+1;
else
    tPow = nPow;
end
[n n2] = subplotn(tPow);

for ipow = 1:nPow
    subplot(n,n2,ipow)
    errorbar([ratios;ratios]', pct_right(:,:,ipow)',(pct_right(:,:,ipow)-squeeze(ci_right(1,:,:,ipow)))',(squeeze(ci_right(2,:,:,ipow))-pct_right(:,:,ipow))', '-o')
    set(gca,'Xscale','log')
    title(['n = ' num2str(sum(totR(1,:,ipow),2)) ' Con, ' num2str(sum(totR(2,:,ipow),2)) ' LED = ' num2str(pows(ipow)) ' mW; ' num2str(length(datelist{ipow})) ' days'])
    hold on
%     text(0.01, .7, num2str(datelist{ipow}'))
    xlabel('Contrast ratio (R/L)')
    ylabel('Fraction Right Choices')
end
if isempty(pow_use) 
    pow_use = pows;
end

if length(pow_use) == nPow 
    tot_all = zeros(2, length(ratios));
    totR_all = zeros(2, length(ratios));
    pct_right_all = zeros(2,length(ratios));
    ci_right_all = zeros(2,2,length(ratios));
    for irat = 1:length(ratios)
        for ib = 1:2
            tot_all(ib,irat) = sum(ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1);
            totR_all(ib,irat) = sum(ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight);
            [pct_right_all(ib,irat) ci_right_all(:,ib,irat)] = binofit(totR_all(ib,irat),tot_all(ib,irat));
        end
    end
    if nPow>1
        subplot(n,n2,tPow)
        errorbar([ratios;ratios]', pct_right_all(:,:)',(pct_right_all(:,:)-squeeze(ci_right_all(1,:,:)))',(squeeze(ci_right_all(2,:,:))-pct_right_all(:,:))', '-o')
        set(gca,'Xscale','log')
        title([mouse '- n = ' num2str(sum(tot_all(1,:),2)) ' Con, ' num2str(sum(tot_all(2,:),2)) ' All LED powers'])
        xlabel('Contrast ratio (R/L)')
        ylabel('Fraction Right Choices')
    end
    suptitle(mouse)
else
    tot_all = zeros(2, length(ratios));
    totR_all = zeros(2, length(ratios));
    pct_right_all = zeros(2,length(ratios));
    ci_right_all = zeros(2,2,length(ratios));
    for irat = 1:length(ratios)
        for ib = 1:2
            for ipow = 1:length(pow_use)
                tot_all(ib,irat) = tot_all(ib,irat) + sum(ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & pow_mat==pow_use(ipow));
                totR_all(ib,irat) = totR_all(ib,irat) + sum(ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight & pow_mat==pow_use(ipow));
                [pct_right_all(ib,irat) ci_right_all(:,ib,irat)] = binofit(totR_all(ib,irat),tot_all(ib,irat));
            end
        end
    end
    if nPow>1
        subplot(n,n2,tPow)
        errorbar([ratios;ratios]', pct_right_all(:,:)',(pct_right_all(:,:)-squeeze(ci_right_all(1,:,:)))',(squeeze(ci_right_all(2,:,:))-pct_right_all(:,:))', '-o')
        set(gca,'Xscale','log')
        title([mouse '- n = ' num2str(sum(tot_all(1,:),2)) ' Con, ' num2str(sum(tot_all(2,:),2)) ': LED powers = ' num2str(pow_use) ' mW'])
        xlabel('Contrast ratio (R/L)')
        ylabel('Fraction Right Choices')
    end
    suptitle(mouse)
end

% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\'...
%     mouse '_Summary.pdf'],'-dpdf','-bestfit'); 

% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\'...
%     mouse '_Summary.pdf'],'-dpdf','-bestfit'); 


tCon_mat = celleqel2mat_padded(input.tGratingContrast);
ind_use =[];
for ipow = 1:length(pow_use)
    ind_use = [ind_use find(pow_mat==pow_use(ipow))];
end
ratio_mat = ratio_mat(ind_use);
b2Ix = b2Ix(ind_use);
IiX = IiX(ind_use);
tRight = tRight(ind_use);
tCon_mat = tCon_mat(ind_use);

% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\'...
%     mouse '_trialInfo.mat'],'ratio_mat', 'b2Ix', 'IiX', 'tRight');

% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\'...
%     mouse '_trialInfo.mat'],'ratio_mat', 'b2Ix', 'IiX', 'tRight');

if doFit
    figure;
    data = zeros(length(ratios),3,2);
    options = struct;
    options.sigmoidName = 'logn';
    options.expType = 'YesNo';
    plotOptions = struct;
    plotOptions.plotData= 0;
    colmat = strvcat('k','g');
    hline = gobjects(1,2);
    for ib = 1:2
        data(:,1,ib) = ratios;
        data(:,2,ib) = totR_all(ib,:)';
        data(:,3,ib) = tot_all(ib,:)';
        result(ib) = psignifit(data(:,:,ib),options);
        if ib == 1
            plotOptions.lineColor= [0 0 0];
        else
            plotOptions.lineColor= [0 1 0];
        end
        [hline(ib),hdata] = plotPsych(result(ib),plotOptions);
        hold on
        errorbar(ratios', pct_right_all(ib,:)',(pct_right_all(ib,:)'-squeeze(ci_right_all(1,ib,:)))',(squeeze(ci_right_all(2,ib,:))...
            -pct_right_all(ib,:)')', ['o' colmat(ib,:)])
    end
%     xlim([0.01 100])
%     axis square
    title([mouse ': ' num2str(pow_use) ' mW'])
    xlabel('Contrast ratio (R/L)')
    ylabel('Fraction Right Choices')
    legend([hline],['Control: n = ' num2str(sum(tot_all(1,:),2))],...
        ['LED: n = ' num2str(sum(tot_all(2,:),2))], 'Location', 'northwest')

%     print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\'...
%         mouse '_Summary_withFit.pdf'],'-dpdf','-bestfit');
%     
%       print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
%         mouse '_Summary_withFit.pdf'],'-dpdf','-bestfit');
%     save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
%         mouse '_fitInfo.mat'],'result', 'ratios', 'totR_all', 'tot_all', 'ib');
end

%% plot each target contrast seperately
% 
% do_tFit = 1;
% 
% edges = [0.0095 0.04 0.2 0.5 1 2.5 5 80 102];
% [bin,n,idx] = histcounts(ratio_mat, edges);
% tRatBin = zeros(size(ratio_mat));
% for i = 1:length(edges)
%     tRatBin(find(idx==i)) = mean(ratio_mat(find(idx==i)));
% end
% 
% t_ratios = unique(tRatBin);
% n_ratios = length(t_ratios);
% tCon = unique(tCon_mat);
% nCon = length(tCon);
% t_totR = zeros(2,length(t_ratios),nCon);
% t_tot_all = zeros(2,length(t_ratios),nCon);
% t_pct_right = zeros(2,length(t_ratios),nCon);
% t_ci_right = zeros(2,2,length(t_ratios),nCon);
% 
% [n,n2] = subplotn(nCon);
% figure;
% colmat = strvcat('k','g');
% for iCon = 1:nCon
%     subplot(n,n2,iCon)
%     for ib = 1:2
%         for irat = 1:n_ratios
%             t_totR(ib,irat,iCon) = sum(tCon_mat==tCon(iCon) & tRatBin==t_ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight);
%             t_tot_all(ib,irat,iCon) = sum(tCon_mat==tCon(iCon) & tRatBin==t_ratios(irat) & IiX==0 & b2Ix==ib-1);
%             [t_pct_right(ib,irat,iCon),t_ci_right(:,ib,irat,iCon)] =...
%                 binofit(sum(tCon_mat==tCon(iCon) & tRatBin==t_ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight),t_tot_all(ib,irat,iCon));
%         end
%         errorbar([t_ratios;t_ratios]', t_pct_right(:,:,iCon)',(t_pct_right(:,:,iCon)-squeeze(t_ci_right(1,:,:,iCon)))',(squeeze(t_ci_right(2,:,:,iCon))-t_pct_right(:,:,iCon))', ['-o' colmat(1)])
%     end
%     set(gca,'Xscale','log')
% %     title(['n = ' num2str(sum(totR(1,:,ipow),2)) ' Con, ' num2str(sum(totR(2,:,ipow),2)) ' LED = ' num2str(pows(ipow)) ' mW; ' num2str(length(datelist{ipow})) ' days'])
%     hold on
%     plot([1,1],[0 1], 'k')
%     xlabel('Contrast ratio (R/L)')
%     ylabel('Fraction Right Choices')
% end
% 
% if do_tFit
%     figure;
%     tdata = zeros(length(t_ratios),3,2,nCon);
%     options = struct;
%     options.sigmoidName = 'logn';
%     options.expType = 'YesNo';
%     plotOptions = struct;
%     plotOptions.plotData= 0;
% %   plotOptions.lineWidth = 1.5;
%     colmat = strvcat('k','g');
%     hline = gobjects(1,2);
%     [n,n2] = subplotn(nCon);
%     for iCon = 1:nCon
%         for ib = 1:2
%             subplot(n,n2,iCon)
%             tdata(:,1,ib,iCon) = t_ratios;
%             tdata(:,2,ib,iCon) = t_totR(ib,:,iCon)';
%             tdata(:,3,ib,iCon) = t_tot_all(ib,:,iCon)';
%             tresult(ib) = psignifit(tdata(:,:,ib,iCon),options);
%             plotOptions.lineColor = colmat(ib);
%             plotOptions.dataColor = colmat(ib);
%             [hline(ib),hdata] = plotPsych(tresult(ib),plotOptions);
%             hold on
%             errorbar(t_ratios', t_pct_right(ib,:,iCon)',(t_pct_right(ib,:,iCon)'-squeeze(t_ci_right(1,ib,:,iCon)))',(squeeze(t_ci_right(2,ib,:,iCon))-t_pct_right(ib,:,iCon)')', ['o' colmat(ib,:)])
%             xlabel('Contrast ratio (R/L)')
%             ylabel('Fraction Right Choices')
%             title(['Target Contrast = ' num2str(chop(tCon(iCon),2))])
%         end
%     end
%     suptitle([mouse ' Sorted by Target Contrast'])
% %     print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\nevyana\Analysis\contrastdiscrim\'...
% %        mouse 'SortedByTargetContrast.pdf'],'-dpdf','-bestfit');
% end

             
% %% fit multiple mice
% 
% mouse_list = strvcat('i459', 'i460', 'i463','i414', 'i461', 'i581');
% nmouse = size(mouse_list,1);
% % data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\nevyana\Analysis\contrastDiscrim\';
% data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\multifit\';
% 
% figure;
% options = struct;
% options.sigmoidName = 'logn';
% options.expType = 'YesNo';
% plotOptions = struct;
% plotOptions.lineWidth = 1.5;
% colmat = char('k','g');
% colmat2 = char('k','b');
% hline = gobjects(1,2);
% tot_control = zeros(1,nmouse);
% tot_LED = zeros(1,nmouse);
% %for imouse = 1:nmouse
% for imouse = 1:3
%     mouse = mouse_list(imouse,:);
%     load(fullfile(data_path,[mouse '_fitInfo.mat']));
%     tot_control(imouse) = sum(tot_all(1,:));
%     tot_LED(imouse) = sum(tot_all(2,:));
%      for ib = 1:2
%         plotOptions.lineColor = colmat(ib);
%         plotOptions.dataColor = colmat(ib);
%         [hline(ib),hdata] = plotPsych(result(ib),plotOptions);
%         hold on
%      end
% end
% 
% for imouse = 4:6
%     mouse = mouse_list(imouse,:);
%     load(fullfile(data_path,[mouse '_fitInfo.mat']));
%     tot_control(imouse) = sum(tot_all(1,:));
%     tot_LED(imouse) = sum(tot_all(2,:));
%      for ib = 1:2
%         plotOptions.lineColor = colmat2(ib);
%         plotOptions.dataColor = colmat2(ib);
%         [hline(ib),hdata] = plotPsych(result(ib),plotOptions);
%         hold on
%      end
% end
% tot_all_control = sum(tot_control);
% tot_all_LED = sum(tot_LED);
% xlabel('Contrast ratio (R/L)')
% ylabel('Fraction Right Choices')
% title(['Contrast Discrimination SOM Arch vs ChR2  n = ' num2str(nmouse) ' mice,' ' Control: n = ' num2str(tot_all_control)])
% legend([hline],['Green LED (Arch): n = ' num2str(sum(tot_LED(1:2)))],...
%         ['Blue LED (ChR2): n = ' num2str(sum(tot_LED(3:4)))],  'Location', 'northwest')
%     
%    % ['Control: n = ' num2str(tot_all_control)],...
% 
% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\multifit\'...
%         'ContrastDiscrim_summary_SOMArchChr2.pdf'],'-dpdf','-bestfit');


%%
% tCon = celleqel2mat_padded(input.tGratingContrast);
% cons = unique(tCon);
% nCon = length(cons);
% totR = zeros(2, length(ratios));
% pct_right = zeros(2,length(ratios));
% ci_right = zeros(2,2,length(ratios));
% figure;
% for icon = 1:nCon
%     subplot(1,4,icon)
%     for irat = 1:length(ratios)
%         for ib = 1:2
%             totR(ib,irat) = sum(tCon==cons(icon) & ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1);
%             [pct_right(ib,irat) ci_right(:,ib,irat)] = binofit(sum(tCon==cons(icon) & ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight),totR(ib,irat));
%         end
%     end
%     errorbar([ratios;ratios]', pct_right',(pct_right-squeeze(ci_right(1,:,:)))',(squeeze(ci_right(2,:,:))-pct_right)', '-o')
%     set(gca,'Xscale','log')
%     title(['Con = ' num2str(chop(cons(icon),2)) '- n = ' num2str(sum(totR(1,:),2)) ' Control, ' num2str(sum(totR(2,:),2)) ' LED'])
%     hold on
%     xlabel('Contrast ratio (R/L)')
%     ylabel('Fraction Right Choices')
% end
% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' mouse '_SummarByTargCon.pdf'],'-dpdf','-bestfit'); 
%             