clear all
close all

%%
behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
mouse = 'i442';
ds= [mouse '_exptList'];
eval(ds)
doFit = 0;
fprintf([mouse '\n'])

%%

clear temp;
good_exp = [];
for id = 1:size(dates,1)
    fprintf([num2str(dates(id,:)) '\n'])
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
        [s b] = selectCalc(input,trials{id});
        if s>=0.9 & abs(b)<0.1
            if ~isempty(trials{id})
                input = trialChopper(input,trials{id});
            end
            if input.doContrastDiscrim
                temp(id).rightGratingContrast = input.rightGratingContrast;
                temp(id).leftGratingContrast = input.leftGratingContrast;
            end
            temp(id).tRightResponse = input.tRightResponse;
            temp(id).tBlock2TrialNumber = input.tBlock2TrialNumber;
            temp(id).trialOutcomeCell = input.trialOutcomeCell;
            temp(id).tGratingContrast = input.tGratingContrast;
            temp(id).tTrialLaserPowerMw = double(input.block2TrialLaserPowerMw).*ones(size(celleqel2mat_padded(input.tTrialLaserPowerMw)));
            temp(id).date = dates(id,:).*ones(size(celleqel2mat_padded(input.tTrialLaserPowerMw)));
            good_exp = [good_exp id];
        end
    end
end

input = concatenateDataBlocks(temp);    
rCon = celleqel2mat_padded(input.rightGratingContrast);
lCon = celleqel2mat_padded(input.leftGratingContrast);
ratio_mat = chop(rCon./lCon,3);
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
            [pct_right(ib,irat,ipow) ci_right(:,ib,irat,ipow)] = binofit(sum(pow_mat==pows(ipow) & ratio_mat==ratios(irat) & IiX==0 & b2Ix==ib-1 & tRight),totR(ib,irat,ipow));
        end
    end
    datelist{ipow} = unique(input.date(find(pow_mat==pows(ipow))));
end
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
    title([mouse '- n = ' num2str(sum(totR(1,:,ipow),2)) ' Control, ' num2str(sum(totR(2,:,ipow),2)) ' LED = ' num2str(pows(ipow)) ' mW'])
    hold on
    text(0.01, .7, num2str(datelist{ipow}'))
    xlabel('Contrast ratio (R/L)')
    ylabel('Fraction Right Choices')
end


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
    title([mouse '- n = ' num2str(sum(tot_all(1,:),2)) ' Control, ' num2str(sum(tot_all(2,:),2)) ' All LED powers'])
    xlabel('Contrast ratio (R/L)')
    ylabel('Fraction Right Choices')
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' mouse '_Summary.pdf'],'-dpdf','-bestfit'); 


if doFit
    figure;
    data = zeros(length(ratios),3,2);
    options = struct;
    options.sigmoidName = 'logn';
    options.expType = 'YesNo';
    plotOptions = struct;
    plotOptions.plotData= 0;
    colmat = strvcat('k','g');
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
        plotPsych(result(ib),plotOptions);
        hold on
        errorbar(ratios', pct_right_all(ib,:)',(pct_right_all(ib,:)'-squeeze(ci_right_all(1,ib,:)))',(squeeze(ci_right_all(2,ib,:))-pct_right_all(ib,:)')', ['o' colmat(ib,:)])
    end
    %xlim([0.01 100])
    axis square
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\2AFC\' mouse '_Summary_withFit.pdf'],'-dpdf','-bestfit'); 
end

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