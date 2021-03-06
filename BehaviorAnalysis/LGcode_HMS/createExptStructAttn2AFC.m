function mouse = createExptStructAttn2AFC;
    rc = behavConstsAttn2AFC;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAttn2AFC;
    
    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).name = pv.mouse_mat(imouse);
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).posN = [];
        mouse(imouse).oddsRightPctN = [];
        for itest = 1:3
            mouse(imouse).test(itest).oddsRightPct = [];
            mouse(imouse).test(itest).testCon = [];
            mouse(imouse).test(itest).intensities = [];
            test(itest).n = 0;
            for ipos = 1:2
                mouse(imouse).test(itest).pos(ipos).nCorrect = [];
                mouse(imouse).test(itest).pos(ipos).nIncorrect = [];
                mouse(imouse).test(itest).pos(ipos).nIgnore = [];
                for ihalf = 1:2
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect = [];
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect = [];
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nIgnore = [];
                end
                for i = 1:2
                    mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimes = cell(1,3);
                    for ihalf = 1:2
                        mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimes = cell(1,3);
                    end
                end
            end
        end
        for iexp = 1:length(mouse(imouse).ind)
            ind = mouse(imouse).ind(iexp);
            fbase = fullfile(rc.pathStr,['data-i' num2str(xd.Subject(ind)) '-' cell2mat(xd.DateStr(ind))]);
            disp(['i' num2str(xd.Subject(ind)) ' ' cell2mat(xd.DateStr(ind))])
            if ~isnan(xd.ChooseMatFile(ind))
                load([fbase '-' num2str(xd.ChooseMatFile(ind)) '.mat'])
            else
                n = dir([fbase '*']);
                if length(n) == 1
                    load(fullfile(rc.pathStr, n.name));
                else
                    error('Need to choose mat file')
                end
            end
            
            trials = str2num(cell2mat(xd.TrialRangeToUse(ind)));
            if isnan(trials)
                trials = input.trialSinceReset;
                range = 1:trials;
            else
                range = trials;
            end
            mouse(imouse).expt(iexp).posN = unique(cell2mat(input.tGratingEccentricityDeg));
            mouse(imouse).expt(iexp).intensities = unique(chop(cell2mat(input.tGratingContrast),3));
            mouse(imouse).expt(iexp).pctCorrect = length(intersect(find(strcmp(input.trialOutcomeCell, 'success')), range))./length(range);
            reactTimes = cell2mat(input.tDecisionTimeMs);
            for ipos = 1:length(mouse(imouse).expt(iexp).posN)
                mouse(imouse).expt(iexp).pos(ipos).trial_ind = intersect(find(cell2mat(input.tGratingEccentricityDeg) == mouse(imouse).expt(iexp).posN(ipos)), range);
                mouse(imouse).expt(iexp).pos(ipos).intensities = unique(chop(unique(cell2mat(input.tGratingContrast(mouse(imouse).expt(iexp).pos(ipos).trial_ind))),3));
            end
            for ipos = 1:2
                for icon = 1:length(mouse(imouse).expt(iexp).intensities)
                    con_ind = find(chop(cell2mat(input.tGratingContrast),3)==mouse(imouse).expt(iexp).intensities(icon));
                    S_ind = intersect(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind), find(strcmp(input.trialOutcomeCell, 'success')));
                    I_ind = intersect(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind), find(strcmp(input.trialOutcomeCell, 'incorrect')));
                    mouse(imouse).expt(iexp).pos(ipos).nTrials(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind));                    
                    mouse(imouse).expt(iexp).pos(ipos).correctInd{icon} = intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'success'))));
                    mouse(imouse).expt(iexp).pos(ipos).incorrectInd{icon} = intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'incorrect'))));
                    mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'success')))));
                    mouse(imouse).expt(iexp).pos(ipos).nIncorrect(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'incorrect')))));
                    mouse(imouse).expt(iexp).pos(ipos).nIgnore(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'ignore')))));
                    mouse(imouse).expt(iexp).pos(ipos).outcome(1).reactTimes{icon} = reactTimes(1,S_ind);
                    mouse(imouse).expt(iexp).pos(ipos).outcome(2).reactTimes{icon} = reactTimes(1,I_ind);
                    mouse(imouse).expt(iexp).pos(ipos).outcome(1).name = 'correct';
                    mouse(imouse).expt(iexp).pos(ipos).outcome(2).name = 'incorrect';
                    [x y] = binofit(mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon),mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon)+mouse(imouse).expt(iexp).pos(ipos).nIncorrect(icon));
                    mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect = x;
                    mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95 = y;
                    for ihalf = 1:2
                        if ihalf == 1
                            sub_range = range(1:ceil(length(range)./2));
                        else
                            sub_range = range(ceil(length(range)./2):end);
                        end
                        S_ind = intersect(sub_range,intersect(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind), find(strcmp(input.trialOutcomeCell, 'success'))));
                        I_ind = intersect(sub_range,intersect(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind), find(strcmp(input.trialOutcomeCell, 'incorrect'))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nTrials(icon) = length(intersect(sub_range, intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind)));                    
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).correctInd{icon} = intersect(sub_range,intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'success')))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).incorrectInd{icon} = intersect(sub_range,intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'incorrect')))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nCorrect(icon) = length(intersect(sub_range,intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'success'))))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nIncorrect(icon) = length(intersect(sub_range,intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'incorrect'))))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nIgnore(icon) = length(intersect(sub_range,intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'ignore'))))));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).outcome(1).reactTimes{icon} = reactTimes(1,S_ind);
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).outcome(2).reactTimes{icon} = reactTimes(1,I_ind);
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).outcome(1).name = 'correct';
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).outcome(2).name = 'incorrect';
                        [x y] = binofit(mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nCorrect(icon),mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nCorrect(icon)+mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nIncorrect(icon));
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).binofit(icon).pctCorrect = x;
                        mouse(imouse).expt(iexp).pos(ipos).half(ihalf).binofit(icon).ci95 = y;
                    end
                end
            end
            if isfield(input, 'doLeftSeparateOdds')
                if input.doLeftSeparateOdds
                    mouse(imouse).expt(iexp).doSepOdds = 1;
                    mouse(imouse).expt(iexp).oddsLeft = [input.leftTrPer80Level1 input.leftTrPer80Level2 input.leftTrPer80Level3 input.leftTrPer80Level4 input.leftTrPer80Level5 input.leftTrPer80Level6 input.leftTrPer80Level7 input.leftTrPer80Level8];   
                    mouse(imouse).expt(iexp).oddsLeftPct = sum(mouse(imouse).expt(iexp).oddsLeft,2)/80;
                    mouse(imouse).expt(iexp).oddsRight = input.trPer80V;
                    mouse(imouse).expt(iexp).oddsRightPct = sum(input.trPer80V,2)/80;
                else
                    mouse(imouse).expt(iexp).doSepOdds = 0;
                    mouse(imouse).expt(iexp).oddsLeft = input.trPer80V/2;
                    mouse(imouse).expt(iexp).oddsRight = input.trPer80V/2;
                    mouse(imouse).expt(iexp).oddsLeftPct = .5;
                    mouse(imouse).expt(iexp).oddsRightPct = .5;
                end
            else
                mouse(imouse).expt(iexp).doSepOdds = 0;
                mouse(imouse).expt(iexp).oddsLeft = input.trPer80V/2;
                mouse(imouse).expt(iexp).oddsRight = input.trPer80V/2;
                mouse(imouse).expt(iexp).oddsLeftPct = .5;
                mouse(imouse).expt(iexp).oddsRightPct = .5;
            end
            if mouse(imouse).expt(iexp).oddsRightPct > mouse(imouse).expt(iexp).oddsLeftPct
                mouse(imouse).expt(iexp).side = 'R';
                itest = 3;
                test_pos = 2;
            elseif mouse(imouse).expt(iexp).oddsRightPct < mouse(imouse).expt(iexp).oddsLeftPct
                mouse(imouse).expt(iexp).side = 'L';
                itest = 1;
                test_pos = 1;
            elseif mouse(imouse).expt(iexp).oddsRightPct == mouse(imouse).expt(iexp).oddsLeftPct
                mouse(imouse).expt(iexp).side = 'B';
                itest = 2;
                test_pos = 1;
            end
            test(itest).n = test(itest).n+1;
            mouse(imouse).test(itest).oddsRightPct = [mouse(imouse).test(itest).oddsRightPct; mouse(imouse).expt(iexp).oddsRightPct];
            mouse(imouse).test(itest).testCon = [mouse(imouse).test(itest).testCon; mouse(imouse).expt(iexp).pos(test_pos).intensities];
            mouse(imouse).test(itest).intensities = [mouse(imouse).test(itest).intensities; mouse(imouse).expt(iexp).intensities];
            mouse(imouse).posN = [mouse(imouse).posN; mouse(imouse).expt(iexp).posN];
            mouse(imouse).oddsRightPctN = [mouse(imouse).oddsRightPctN, mouse(imouse).expt(iexp).oddsRightPct];
            for ipos = 1:2
                mouse(imouse).test(itest).pos(ipos).nCorrect = [mouse(imouse).test(itest).pos(ipos).nCorrect; mouse(imouse).expt(iexp).pos(ipos).nCorrect];
                mouse(imouse).test(itest).pos(ipos).nIncorrect = [mouse(imouse).test(itest).pos(ipos).nIncorrect; mouse(imouse).expt(iexp).pos(ipos).nIncorrect];
                mouse(imouse).test(itest).pos(ipos).nIgnore = [mouse(imouse).test(itest).pos(ipos).nIgnore; mouse(imouse).expt(iexp).pos(ipos).nIgnore];
                for ihalf = 1:2
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect = [mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect; mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nCorrect];
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect = [mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect; mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nIncorrect];
                    mouse(imouse).test(itest).pos(ipos).half(ihalf).nIgnore = [mouse(imouse).test(itest).pos(ipos).half(ihalf).nIgnore; mouse(imouse).expt(iexp).pos(ipos).half(ihalf).nIgnore];
                end
                for icon = 1:length(mouse(imouse).expt(iexp).intensities)
                    for i = 1:2
                        mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimes{test(itest).n,icon} = mouse(imouse).expt(iexp).pos(ipos).outcome(i).reactTimes{icon};
                        for ihalf = 1:2
                        	mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimes{test(itest).n,icon} = mouse(imouse).expt(iexp).pos(ipos).half(ihalf).outcome(i).reactTimes{icon};
                        end
                    end
                end
                    mouse(imouse).test(itest).pos(ipos).outcome(1).name = 'correct';
                    mouse(imouse).test(itest).pos(ipos).outcome(2).name = 'incorrect';
            end
        end
        intensities_all = [];
        for itest = 1:3
            intensities_all = [intensities_all; mouse(imouse).test(itest).intensities];
        end
        mouse(imouse).intensities = unique(intensities_all);
        for itest = 1:3
            for ipos = 1:2
                if length(mouse(imouse).test(itest).pos(ipos).nCorrect + mouse(imouse).test(itest).pos(ipos).nIncorrect)>1
                    for icon = 1:length(mouse(imouse).intensities)
                        ind = find(mouse(imouse).test(itest).intensities == mouse(imouse).intensities(icon));
                        if length(ind)>1
                            [x y] = binofit(sum(mouse(imouse).test(itest).pos(ipos).nCorrect(ind),1),sum(mouse(imouse).test(itest).pos(ipos).nCorrect(ind),1) + sum(mouse(imouse).test(itest).pos(ipos).nIncorrect(ind),1));
                            mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect = x;
                            mouse(imouse).test(itest).pos(ipos).con(icon).ci95 = y;
                            for i = 1:2
                                mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimesByCon{icon} = [];
                                for iind = 1:length(ind)
                                    ii = ind(iind);
                                    mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimesByCon{icon} = [mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimesByCon{icon} mouse(imouse).test(itest).pos(ipos).outcome(i).reactTimes{ii}];
                                end
                            end
                            for ihalf = 1:2
                                [x y] = binofit(sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect(ind),1),sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nCorrect(ind),1) + sum(mouse(imouse).test(itest).pos(ipos).half(ihalf).nIncorrect(ind),1));
                                mouse(imouse).test(itest).pos(ipos).half(ihalf).con(icon).pctCorrect = x;
                                mouse(imouse).test(itest).pos(ipos).half(ihalf).con(icon).ci95 = y;
                                for i = 1:2
                                    mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimesByCon{icon} = [];
                                    for iind = 1:length(ind)
                                        ii = ind(iind);
                                        mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimesByCon{icon} = [mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimesByCon{icon} mouse(imouse).test(itest).pos(ipos).half(ihalf).outcome(i).reactTimes{ii}];
                                    end
                                end
                            end
                        else
                             mouse(imouse).test(itest).pos(ipos).con(icon).pctCorrect = [];
                             for ihalf = 1:2
                                 mouse(imouse).test(itest).pos(ipos).half(ihalf).con(icon).pctCorrect = [];
                             end
                        end
                    end
                end
            end
        end
    end
    save([rc.fitOutputSummary '\' date '-2afc-attn_mouse.mat'], 'mouse')
end

    
    