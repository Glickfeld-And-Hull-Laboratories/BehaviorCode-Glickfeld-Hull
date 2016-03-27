rc = behavConstsAttn2AFC;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
pv = behavParamsAttn2AFC;

colmat = strvcat('g', 'k', 'b');

%% avg hit rate by contrast
c_all = cell(1, size(pv.mouse_mat,2));
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    ntest = zeros(3,3,2);
    c = zeros(3,length(mouse(imouse).intensities));
    for itest = 1:length(tests)
        for icon = 1:length(mouse(imouse).intensities)
            c(itest,icon) = ~isempty(mouse(imouse).test(itest).pos(1).con(icon).pctCorrect);
        end
    end
    [n, n2] = subplotn(max(sum(c,2),[],1));
    indn = [];
    c_all{imouse} = c;
    for itest = 1:length(tests)
        num(itest) = length(mouse(imouse).test(itest).oddsRightPct);
        for icon = 1:length(mouse(imouse).intensities)
        	if c(itest,icon)
                indc = find(sum(c(:,1:icon),1)>0);
                indn = [indn icon];
                subplot(n,n2,length(indc))
                x = mouse(imouse).posN(1,:);
                y = [mouse(imouse).test(itest).pos(1).con(icon).pctCorrect mouse(imouse).test(itest).pos(2).con(icon).pctCorrect];
                errL = [mouse(imouse).test(itest).pos(1).con(icon).pctCorrect - mouse(imouse).test(itest).pos(1).con(icon).ci95(:,1) mouse(imouse).test(itest).pos(2).con(icon).pctCorrect - mouse(imouse).test(itest).pos(2).con(icon).ci95(:,1)];
                errU = [mouse(imouse).test(itest).pos(1).con(icon).ci95(:,2) - mouse(imouse).test(itest).pos(1).con(icon).pctCorrect mouse(imouse).test(itest).pos(2).con(icon).ci95(:,2) - mouse(imouse).test(itest).pos(2).con(icon).pctCorrect];
                errorbar(x, y, errL, errU, ['-o' colmat(itest,:)]);
                hold on
                ind = find(mouse(imouse).test(itest).intensities == mouse(imouse).intensities(icon));
                for ipos = 1:2
                    ntest(icon,itest,ipos) = sum(mouse(imouse).test(itest).pos(ipos).nCorrect(ind))+sum(mouse(imouse).test(itest).pos(ipos).nIncorrect(ind));
                end
            end
        end
    end
    uindn = unique(indn);
    for icon = 1:length(uindn)
        con = uindn(icon);
        subplot(n,n2,icon) 
        title({[num2str(mouse(imouse).intensities(con)*100) '% contrast'], ['Black (' num2str(ntest(con,2,:)) '); Green (' num2str(ntest(con,1,:)) '); Blue (' num2str(ntest(con,3,:)) ')']})
        ylim([0 1])
        ylabel('Hit rate')
        xlabel('Stimulus position (deg)')
    end
    %legend([num2str((tests*100)') repmat('% 30 deg- n = ', 3,1) num2str(num')])
    suptitle([num2str(pv.mouse_mat(:,imouse)) '- Black- 50% - n = ' num2str(num(2)) '; Green- 10% - n = ' num2str(num(1)) '; Blue- 90% - n = ' num2str(num(1))])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_avgHitRate_probeCon.pdf'], '-dpdf')
end

%% avg hit rate by contrast by half
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    c = c_all{imouse};
    tests = unique(mouse(imouse).oddsRightPctN);
    ntest = zeros(length(mouse(imouse).intensities),3,2,2);
    n = (max(sum(c,2),[],1));
    indn = [];
    ncon = length(find(sum(c,1)>0));
    for itest = 1:length(tests)
        num(itest) = length(mouse(imouse).test(itest).oddsRightPct);
        for icon = 1:length(mouse(imouse).intensities)
        	if c(itest,icon)
                for ihalf = 1:2
                    indc = find(sum(c(:,1:icon),1)>0);
                    indn = [indn icon];
                    subplot(2,n,length(indc)+((ihalf-1)*ncon))
                    x = mouse(imouse).posN(1,:);
                    y = [mouse(imouse).test(itest).pos(1).half(ihalf).con(icon).pctCorrect mouse(imouse).test(itest).pos(2).half(ihalf).con(icon).pctCorrect];
                    errL = [mouse(imouse).test(itest).pos(1).half(ihalf).con(icon).pctCorrect - mouse(imouse).test(itest).pos(1).half(ihalf).con(icon).ci95(:,1) mouse(imouse).test(itest).pos(2).half(ihalf).con(icon).pctCorrect - mouse(imouse).test(itest).pos(2).half(ihalf).con(icon).ci95(:,1)];
                    errU = [mouse(imouse).test(itest).pos(1).half(ihalf).con(icon).ci95(:,2) - mouse(imouse).test(itest).pos(1).half(ihalf).con(icon).pctCorrect mouse(imouse).test(itest).pos(2).half(ihalf).con(icon).ci95(:,2) - mouse(imouse).test(itest).pos(2).half(ihalf).con(icon).pctCorrect];
                    errorbar(x, y, errL, errU, ['-o' colmat(itest,:)]);
                    hold on
                    ind = find(mouse(imouse).test(itest).intensities == mouse(imouse).intensities(icon));
                    for ipos = 1:2
                        ntest(icon,itest,ipos,ihalf) = sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect(ind),1)+sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect(ind),1);
                    end
                end
            end
        end
    end
    uindn = unique(indn);
    for icon = 1:length(uindn)
        for ihalf = 1:2
            con = uindn(icon);
            subplot(2,n,icon+((ihalf-1)*ncon)) 
            title({[num2str(mouse(imouse).intensities(con)*100) '% contrast'], ['Black (' num2str(ntest(con,2,:,ihalf)) '); Green (' num2str(ntest(con,1,:,ihalf)) '); Blue (' num2str(ntest(con,3,:,ihalf)) ')']})
            ylim([0 1])
            ylabel('Hit rate')
            xlabel('Stimulus position (deg)')
        end
    end
    %legend([num2str((tests*100)') repmat('% 30 deg- n = ', 3,1) num2str(num')])
    suptitle([num2str(pv.mouse_mat(:,imouse)) '- Black- 50% - n = ' num2str(num(2)) '; Green- 10% - n = ' num2str(num(1)) '; Blue- 90% - n = ' num2str(num(1))])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_avgHitRate_probeCon_byHalf.pdf'], '-dpdf')
end


%% psych curves for each day
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
        
%% summary psych curves across days
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    ncon = length(mouse(imouse).intensities);
    intensities = reshape([fliplr(mouse(imouse).intensities)*-1 mouse(imouse).intensities], ncon*2, 1);
    for itest = 1:length(tests)
        subplot(3,1,itest)
        pctCorr = zeros(ncon, 2);
        ci95 = zeros(ncon,2,2);
        for icon = 1:ncon
            for ipos = 1:2
                if isempty(mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect)
                    mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect = NaN;
                    mouse(imouse).test(itest).pos(ipos).con(icon).ci95 = [NaN NaN];
                end
                pctCorr(icon,ipos) = mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect;
                ci95(icon,:,ipos) = mouse(imouse).test(itest).pos(ipos).con(icon).ci95;
            end
        end
        pctCorr(:,1) = 1- pctCorr(:,1);
        ci95(:,:,1) = fliplr(1-ci95(:,:,1));
        pctCorr = reshape(pctCorr, ncon*2,1);
        ci95 = reshape(permute(ci95,[1 3 2]), ncon*2,2);
        errorbar(intensities, pctCorr, pctCorr-ci95(:,1), ci95(:,2)-pctCorr, ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Fract. right choice')
        ylim([0 1])
        xlim([-1.1 1.1])
        title([num2str(tests(itest)*100) ' % right trials'])
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Summary Fraction Right Choice'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_fractRightSummary.pdf'], '-dpdf')
end

figure;
for imouse = 1:size(pv.mouse_mat,2)
    [n, n2] = subplotn(size(pv.mouse_mat,2));
    subplot(n, n2, imouse)
    tests = unique(mouse(imouse).oddsRightPctN);
    ncon = length(mouse(imouse).intensities);
    intensities = reshape([fliplr(mouse(imouse).intensities)*-1 mouse(imouse).intensities], ncon*2, 1);
    for itest = 1:length(tests)
        pctCorr = zeros(ncon, 2);
        ci95 = zeros(ncon,2,2);
        for icon = 1:ncon
            for ipos = 1:2
                if isempty(mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect)
                    mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect = NaN;
                    mouse(imouse).test(itest).pos(ipos).con(icon).ci95 = [NaN NaN];
                end
                pctCorr(icon,ipos) = mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect;
                ci95(icon,:,ipos) = mouse(imouse).test(itest).pos(ipos).con(icon).ci95;
            end
        end
        pctCorr(:,1) = 1- pctCorr(:,1);
        ci95(:,:,1) = fliplr(1-ci95(:,:,1));
        pctCorr = reshape(pctCorr, ncon*2,1);
        ci95 = reshape(permute(ci95,[1 3 2]), ncon*2,2);
        errorbar(intensities, pctCorr, pctCorr-ci95(:,1), ci95(:,2)-pctCorr, ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Fract. right choice')
        ylim([0 1])
        xlim([-1.1 1.1])
        title(['i' num2str(pv.mouse_mat(:,imouse))])
    end
end
print([rc.fitOutputSummary '\' date '_allMice_fractRightSummary.pdf'], '-dpdf')


%% number of ignored trials by contrast across days
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    c = c_all{imouse};
    for itest = 1:length(tests)
        subplot(3,1,itest)
        for icon = 1:length(mouse(imouse).intensities)
            ind = find(mouse(imouse).test(itest).intensities == mouse(imouse).intensities(icon));
            if length(ind)>0
                if c(itest, icon)
                    mouse(imouse).test(itest).ExpNIgnore(icon,:) = [sum(mouse(imouse).test(itest).pos(1).nIgnore(ind),1) sum(mouse(imouse).test(itest).pos(2).nIgnore(ind),1)];
                else
                    mouse(imouse).test(itest).ExpNIgnore(icon,:) = [NaN NaN];
                end
            else
                mouse(imouse).test(itest).ExpNIgnore(icon,:) = [NaN NaN];
            end
        end
        plot([mouse(imouse).intensities*-1; mouse(imouse).intensities], [mouse(imouse).test(itest).ExpNIgnore(:,1); mouse(imouse).test(itest).ExpNIgnore(:,2)] , ['o' colmat(itest,:)]);
        hold on
        xlabel('R-L contrast')
        ylabel('Ignored trials')
        xlim([-1.1 1.1])
        title([num2str(tests(itest)*100) ' % right trials'])
    end
    suptitle([num2str(pv.mouse_mat(:,imouse)) ' Summary Ignores by Target Contrast'])
    print([rc.fitOutputSummary '\' date '_' num2str(pv.mouse_mat(:,imouse)) '_IgnoreSummary.pdf'], '-dpdf')
end

%% react times by contrast and outcome
for imouse = 1:size(pv.mouse_mat,2)
    figure;
    tests = unique(mouse(imouse).oddsRightPctN);
    c = c_all{imouse};
    start = 1;
    for itest = 1:length(tests)
        for i = 1:3
            subplot(3,3,start)
            for ipos = 1:2
                for icon = 1:length(mouse(imouse).intensities)
                    ind = find(mouse(imouse).test(itest).intensities == mouse(imouse).intensities(icon));
                    if length(ind)>0
                        if c(itest, icon)
                            if i < 3
                                reactTimes = mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimesByCon{icon};
                            else
                                reactTimes = [mouse(imouse).test(itest).pos(ipos).outcome(1).reactTimesByCon{icon} mouse(imouse).test(itest).pos(ipos).outcome(2).reactTimesByCon{icon}];
                            end
                            if ipos == 1
                                x = mouse(imouse).intensities(icon)*-1; 
                            else
                                x = mouse(imouse).intensities(icon);
                            end
                            errorbar(x, mean(reactTimes,2), std(reactTimes,[],2)./sqrt(size(reactTimes,2)), ['o' colmat(itest,:)]);
                            hold on
                        end
                    end
                end
            end
            start = start+1;
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
                     