rc = behavConstsAttn2AFC;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
pv = behavParamsAttn2AFC;

colmat = strvcat('g', 'k', 'b');

%avg hit rate by contrast
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    ntest = zeros(3,3,2);
    for itest = 1:length(tests)
        num(itest) = length(mouse(imouse).test(itest).oddsRightPct);
        if size(mouse(imouse).test(itest).pos(1).pctCorrect,1) == 1
            for icon = 1:length(mouse(imouse).test(itest).pos(1).pctCorrect)
                subplot(2,2,icon)
                x = mouse(imouse).posN(1,:);
                y = [mouse(imouse).test(itest).pos(1).pctCorrect(icon) mouse(imouse).test(itest).pos(2).pctCorrect(icon)];
                errL = [mouse(imouse).test(itest).pos(1).pctCorrect(icon) - mouse(imouse).test(itest).pos(1).ci95(icon,1) mouse(imouse).test(itest).pos(2).pctCorrect(icon) - mouse(imouse).test(itest).pos(2).ci95(icon,1)];
                errU = [mouse(imouse).test(itest).pos(1).ci95(icon,2) - mouse(imouse).test(itest).pos(1).pctCorrect(icon) mouse(imouse).test(itest).pos(2).ci95(icon,2) - mouse(imouse).test(itest).pos(2).pctCorrect(icon)];
                errorbar(x, y, errL, errU, ['-o' colmat(itest,:)]);
                hold on
                for ipos = 1:2
                    ntest(:,itest,ipos) = sum(mouse(imouse).test(itest).pos(ipos).nCorrect,1)+sum(mouse(imouse).test(itest).pos(ipos).nIncorrect,1);
                end
                if itest == 3
                    title({[num2str(mouse(imouse).test(itest).intensities(1,icon)*100) '% contrast'], ['Black (' num2str(ntest(icon,2,:)) '); Green (' num2str(ntest(icon,1,:)) '); Blue (' num2str(ntest(icon,3,:)) ')']})
                    ylim([0 1])
                    ylabel('Hit rate')
                    xlabel('Stimulus position (deg)')
                end
            end
        end
    end
    %legend([num2str((tests*100)') repmat('% 30 deg- n = ', 3,1) num2str(num')])
    suptitle([num2str(pv.mouse_mat(:,imouse)) '- Black- 50% - n = ' num2str(num(2)) '; Green- 10% - n = ' num2str(num(1)) '; Blue- 90% - n = ' num2str(num(1))])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_avgHitRate_probeCon.pdf'], '-dpdf')
end

%avg hit rate by contrast by half
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    ntest = zeros(3,3,2);
    for itest = 1:length(tests)
        num(itest) = length(mouse(imouse).test(itest).oddsRightPct);
        if size(mouse(imouse).test(itest).pos(1).half(1).pctCorrect,1) == 1
            ncon = length(mouse(imouse).test(itest).pos(1).pctCorrect);
            for icon = 1:length(mouse(imouse).test(itest).pos(1).pctCorrect)
                for ihalf = 1:2
                    subplot(2,3,icon+((ihalf-1)*ncon))
                    x = mouse(imouse).posN(1,:);
                    y = [mouse(imouse).test(itest).pos(1).half(ihalf).pctCorrect(icon) mouse(imouse).test(itest).pos(2).half(ihalf).pctCorrect(icon)];
                    errL = [mouse(imouse).test(itest).pos(1).half(ihalf).pctCorrect(icon) - mouse(imouse).test(itest).pos(1).half(ihalf).ci95(icon,1) mouse(imouse).test(itest).pos(2).half(ihalf).pctCorrect(icon) - mouse(imouse).test(itest).pos(2).half(ihalf).ci95(icon,1)];
                    errU = [mouse(imouse).test(itest).pos(1).half(ihalf).ci95(icon,2) - mouse(imouse).test(itest).pos(1).half(ihalf).pctCorrect(icon) mouse(imouse).test(itest).pos(2).half(ihalf).ci95(icon,2) - mouse(imouse).test(itest).pos(2).half(ihalf).pctCorrect(icon)];
                    errorbar(x, y, errL, errU, ['-o' colmat(itest,:)]);
                    hold on
                    for ipos = 1:2
                        ntest(:,itest,ipos) = sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect,1)+sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect,1);
                    end
                    if itest == 3
                        title({[num2str(mouse(imouse).test(itest).intensities(1,icon)*100) '% contrast'], ['Black (' num2str(ntest(icon,2,:)) '); Green (' num2str(ntest(icon,1,:)) '); Blue (' num2str(ntest(icon,3,:)) ')']})
                        ylim([0 1])
                        ylabel('Hit rate')
                        xlabel('Stimulus position (deg)')
                    end
                end
            end
        end
    end
    %legend([num2str((tests*100)') repmat('% 30 deg- n = ', 3,1) num2str(num')])
    suptitle([num2str(pv.mouse_mat(:,imouse)) '- Black- 50% - n = ' num2str(num(2)) '; Green- 10% - n = ' num2str(num(1)) '; Blue- 90% - n = ' num2str(num(1))])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_avgHitRate_probeCon_byHalf.pdf'], '-dpdf')
end


%psych curves for each day
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    nexp = length(mouse(imouse).ind);
    n = ceil(sqrt(nexp));
    if (n^2)-n > nexp
        n2 = n-1;
    else
        n2 = n;
    end
    for iexp = 1:nexp
        subplot(n, n2, iexp)
        for ipos = 1:2
            for icon = 1:length(mouse(imouse).expt(iexp).intensities)
                if ipos == 1
                    errorbar(mouse(imouse).expt(iexp).intensities(icon).*(-1), 1- mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect, mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect- mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95(2), mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95(1) - mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect, 'ok')
                    hold on
                else
                    errorbar(mouse(imouse).expt(iexp).intensities(icon), mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect, mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect- mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95(1), mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95(2) - mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect, 'ok')
                    hold on
                end
            end
        end
        ylim([0 1])
        xlim([-1.1 1.1])
        title([num2str(mouse(imouse).expt(iexp).oddsRightPct * 100) ' % Right trials'])
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Fraction Right Choice'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_fractRight_byday.pdf'], '-dpdf')
end

%timecourse of %correct for 12.5% contrast for each day
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    nexp = length(mouse(imouse).ind);
    n = ceil(sqrt(nexp));
    if (n^2)-n > nexp
        n2 = n-1;
    else
        n2 = n;
    end
    for iexp = 1:nexp
        subplot(n, n2, iexp)
        for ipos = 1:2
            for icon = 2
                plotData = NaN(1,max([mouse(imouse).expt(iexp).pos(ipos).correctInd{icon} mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon}],[],2));
                plotData(mouse(imouse).expt(iexp).pos(ipos).correctInd{icon}) = 1;
                plotData(mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon}) = 0;
                smPlotData = smooth(plotData,20);
                mouse(imouse).pos(ipos).plotData_all{iexp} = plotData;
                mouse(imouse).pos(ipos).ntrials(iexp) = length(plotData);
                if ipos ==1
                    plot(mouse(imouse).expt(iexp).pos(ipos).correctInd{icon}, ones(size(mouse(imouse).expt(iexp).pos(ipos).correctInd{icon})), 'xg')
                    hold on
                    plot(mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon}, zeros(size(mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon})), 'xg')
                    hold on
                    plot(1:length(plotData), smPlotData, '-g')
                elseif ipos ==2
                    plot(mouse(imouse).expt(iexp).pos(ipos).correctInd{icon}, ones(size(mouse(imouse).expt(iexp).pos(ipos).correctInd{icon})), 'xb')
                    hold on
                    plot(mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon}, zeros(size(mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon})), 'xb')
                    hold on
                    plot(1:length(plotData), smPlotData, '-b')
                end
            end
        end
        ylim([-1 2])
        title([num2str(mouse(imouse).expt(iexp).oddsRightPct * 100) ' % Right trials'])
        ylabel('1 = Correct')
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' TC corrects'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_TCcorrect_byday.pdf'], '-dpdf')
end

for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    for itest = 1:length(tests)
        ind = find(mouse(imouse).oddsRightPctN == tests(:,itest));
        for ipos = 1:2
            max_trials = max(mouse(imouse).pos(ipos).ntrials(ind),[],2);
            plotData_all = NaN(max_trials, length(ind));
            for i = 1:length(ind)
                plotData_all(1:mouse(imouse).pos(ipos).ntrials(ind(i)),i) = mouse(imouse).pos(ipos).plotData_all{ind(i)};
            end
            plotData_avg = nanmean(plotData_all,2);
            isnan_plotdata = ~isnan(plotData_all);
            plotData_sem = std(plotData_all,[],2)./(sqrt(sum(isnan_plotdata,2)));
            subplot(1,3,itest)
            if ipos == 1;
                errorbar(1:max_trials, plotData_avg', plotData_sem', 'g')
                hold on
            elseif ipos == 2;
                errorbar(1:max_trials, plotData_avg', plotData_sem', 'b')
                hold on
            end
        end
        xlim([0 300])
        ylim([0 1])
        xlabel('Trial number')
        ylabel('Percent correct')
        title([num2str(mouse(imouse).expt(iexp).oddsRightPct * 100) ' % Right trials'])
    end
end

%summary percent correct by test
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    for itest = 1:length(tests)
        ind = find(mouse(imouse).oddsRightPctN == tests(:,itest));
        pctCorrect_all = [];
        for i = 1:length(ind)
            pctCorrect_all = [pctCorrect_all mouse(imouse).expt(i).pctCorrect];
        end
        pctCorrect_avg = mean(pctCorrect_all,2);
        pctCorrect_sem = std(pctCorrect_all,[],2)./sqrt(size(pctCorrect_all,2));
        errorbar(itest,pctCorrect_avg,pctCorrect_sem,['o' colmat(itest)]);
        hold on
    end
    ylim([0 1])
    ylabel('% Correct')
    title('% Correct by condition')
end
        

%summary psych curves across days
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    for itest = 1:length(tests)
        ind = find(mouse(imouse).oddsRightPctN == tests(:,itest));
        subplot(3,1,itest)
        mouse(imouse).test(itest).ExpIntensities = [];
        mouse(imouse).test(itest).ExpNCorrect = 0;
        mouse(imouse).test(itest).ExpNIncorrect = 0;
        for iexp = 1:length(ind)
            i = ind(iexp);
            mouse(imouse).test(itest).ExpIntensities = [mouse(imouse).test(itest).ExpIntensities ; [fliplr(mouse(imouse).expt(i).intensities)*-1 mouse(imouse).expt(i).intensities]];
            mouse(imouse).test(itest).ExpNCorrect = [mouse(imouse).test(itest).ExpNCorrect + [fliplr(mouse(imouse).expt(i).pos(1).nCorrect) mouse(imouse).expt(i).pos(2).nCorrect]];
            mouse(imouse).test(itest).ExpNIncorrect = [mouse(imouse).test(itest).ExpNIncorrect + [fliplr(mouse(imouse).expt(i).pos(1).nIncorrect) mouse(imouse).expt(i).pos(2).nIncorrect]];
        end
        [mouse(imouse).test(itest).pctCorr mouse(imouse).test(itest).ci95] = binofit(mouse(imouse).test(itest).ExpNCorrect,mouse(imouse).test(itest).ExpNCorrect + mouse(imouse).test(itest).ExpNIncorrect);
        xval = mean(mouse(imouse).test(itest).ExpIntensities,1);
        mouse(imouse).test(itest).pctCorr(:,1:3) = 1-mouse(imouse).test(itest).pctCorr(:,1:3);
        mouse(imouse).test(itest).ci95(1:3,:) = fliplr(1-mouse(imouse).test(itest).ci95(1:3,:));
        errorbar(xval, mouse(imouse).test(itest).pctCorr, mouse(imouse).test(itest).pctCorr-mouse(imouse).test(itest).ci95(:,1)', mouse(imouse).test(itest).ci95(:,2)' - mouse(imouse).test(itest).pctCorr, ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Fract. right choice')
        ylim([0 1])
        xlim([-1.1 1.1])
        title([num2str(tests(itest)*100) ' % right trials - n = ' num2str(length(ind)) 'sessions'])
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Summary Fraction Right Choice'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_fractRightSummary.pdf'], '-dpdf')
end

%number of ignored trials by contrast across days
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    for itest = 1:length(tests)
        subplot(3,1,itest)
        mouse(imouse).test(itest).ExpNIgnore = [fliplr(sum(mouse(imouse).test(itest).pos(1).nIgnore,1)) sum(mouse(imouse).test(itest).pos(2).nIgnore,1)];
        xval = mouse(imouse).test(itest).ExpIntensities(1,:);
        plot(xval, mouse(imouse).test(itest).ExpNIgnore, ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Ignored trials')
        xlim([-1.1 1.1])
        title([num2str(tests(itest)*100) ' % right trials - n = ' num2str(length(ind)) 'sessions'])
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Summary Ignores by Target Contrast'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_IgnoreSummary.pdf'], '-dpdf')
end

%react times by contrast and outcome
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    for itest = 1:length(tests)
        xval = mouse(imouse).test(itest).ExpIntensities(1,:);
        for i = 1:3
            subplot(3,3,3*(itest-1)+i)
            start = 1;
            for ipos = 1:2
                if ipos == 1
                    cons = [3 2 1];
                else
                    cons = 1:3;
                end
                for ii = 1:length(cons)
                    icon = cons(ii);
                    if i<3
                        reactTimes = mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimes{icon};
                    else
                        reactTimes = [mouse(imouse).test(itest).pos(ipos).outcome(1).reactTimes{icon} mouse(imouse).test(itest).pos(ipos).outcome(1).reactTimes{icon}];
                    end
                    if length(reactTimes)<1
                        reactTimes = NaN;
                    end
                    errorbar(xval(start), mean(reactTimes,2), std(reactTimes,[],2)./sqrt(size(reactTimes,2)), ['o' colmat(itest,:)]);
                    hold on
                    start = start+1;
                end
            end
            xlabel('R-L contrast')
            ylabel('React times (ms)')
            xlim([-1.1 1.1])
            ylim([0 7000])
            if i == 1
                str_out = 'correct';
            elseif i == 2
                str_out = 'incorrect';
            elseif i == 3
                str_out = 'all';
            end
            title([num2str(tests(itest)*100) ' % right - ' str_out ' trials'])
        end
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Summary React Times by Target Contrast'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_ReactTimeSummary.pdf'], '-dpdf')
end

for imouse = 1 :size(pv.mouse_mat,2)
    itest = 2;
    figure;
    ind2 = find(mouse(imouse).oddsRightPctN == tests(:,itest));
    mouse(imouse).test(itest).intensities = [];
    mouse(imouse).test(itest).nCorrect = 0;
    mouse(imouse).test(itest).nIncorrect = 0;
    for ipos = 1:2
        for iexp = 1:length(ind2)
            i = ind2(iexp);
            mouse(imouse).test(itest).intensities = [mouse(imouse).test(itest).intensities ; mouse(imouse).expt(i).intensities];
            mouse(imouse).test(itest).nCorrect = [mouse(imouse).test(itest).nCorrect + mouse(imouse).expt(i).pos(ipos).nCorrect];
            mouse(imouse).test(itest).nIncorrect = [mouse(imouse).test(itest).nIncorrect + mouse(imouse).expt(i).pos(ipos).nIncorrect];
        end
        [mouse(imouse).test(itest).pctCorr mouse(imouse).test(itest).ci95] = binofit(mouse(imouse).test(itest).nCorrect,mouse(imouse).test(itest).nCorrect + mouse(imouse).test(itest).nIncorrect);
        xval = mean(mouse(imouse).test(itest).intensities,1);
        if ipos == 1
            mouse(imouse).test(itest).pctCorr = 1-mouse(imouse).test(itest).pctCorr;
            mouse(imouse).test(itest).ci95 = fliplr(1-mouse(imouse).test(itest).ci95);
        end
        subplot(1,2,ipos)
        errorbar(xval, mouse(imouse).test(itest).pctCorr, mouse(imouse).test(itest).pctCorr-mouse(imouse).test(itest).ci95(:,1)', mouse(imouse).test(itest).ci95(:,2)' - mouse(imouse).test(itest).pctCorr, ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Fract. right choice')
        ylim([0 1])
        xlim([0.01 2])
        title([num2str(mouse(imouse).posN(1,ipos)) 'deg - n = ' num2str(length(ind2)) 'sessions'])
        set(gca, 'xscale', 'log')
    end
end

            

            