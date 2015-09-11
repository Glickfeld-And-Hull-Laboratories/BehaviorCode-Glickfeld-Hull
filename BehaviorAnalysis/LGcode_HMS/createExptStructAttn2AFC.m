function mouse = createExptStructAttn2AFC;
    rc = behavConstsAttn2AFC;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
    pv = behavParamsAttn2AFC;
    
    for imouse = 1:length(pv.mouse_mat)
        mouse(imouse).ind = find(xd.Subject ==  pv.mouse_mat(imouse));
        mouse(imouse).posN = [];
        mouse(imouse).oddsRightPctN = [];
        for itest = 1:3
            mouse(imouse).test(itest).oddsRightPct = [];
            mouse(imouse).test(itest).testCon = [];
            mouse(imouse).test(itest).intensities = [];
            for ipos = 1:2
                mouse(imouse).test(itest).pos(ipos).nCorrect = [];
                mouse(imouse).test(itest).pos(ipos).nIncorrect = [];
            end
        end
        for iexp = 1:length(mouse(imouse).ind)
            ind = mouse(imouse).ind(iexp);
            fbase = fullfile(rc.pathStr,['data-i' num2str(xd.Subject(ind)) '-' cell2mat(xd.DateStr(ind))]);
            if ~isnan(xd.ChooseMatFile(ind))
                load([fbase '-' num2str(xd.ChooseMatFile) '.mat'])
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
            for ipos = 1:length(mouse(imouse).expt(iexp).posN)
                mouse(imouse).expt(iexp).pos(ipos).trial_ind = intersect(find(cell2mat(input.tGratingEccentricityDeg) == mouse(imouse).expt(iexp).posN(ipos)), range);
                mouse(imouse).expt(iexp).pos(ipos).intensities = unique(chop(unique(cell2mat(input.tGratingContrast(mouse(imouse).expt(iexp).pos(ipos).trial_ind))),3));
            end
            for ipos = 1:2
                for icon = 1:length(mouse(imouse).expt(iexp).intensities)
                    con_ind = find(chop(cell2mat(input.tGratingContrast),3)==mouse(imouse).expt(iexp).intensities(icon));
                    mouse(imouse).expt(iexp).pos(ipos).nTrials(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, con_ind));                    
                    mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'success')))));
                    mouse(imouse).expt(iexp).pos(ipos).nIncorrect(icon) = length(intersect(mouse(imouse).expt(iexp).pos(ipos).trial_ind, intersect(con_ind, find(strcmp(input.trialOutcomeCell, 'incorrect')))));
                    [x y] = binofit(mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon),mouse(imouse).expt(iexp).pos(ipos).nCorrect(icon)+mouse(imouse).expt(iexp).pos(ipos).nIncorrect(icon));
                    mouse(imouse).expt(iexp).pos(ipos).binofit(icon).pctCorrect = x;
                    mouse(imouse).expt(iexp).pos(ipos).binofit(icon).ci95 = y;
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
            mouse(imouse).test(itest).oddsRightPct = [mouse(imouse).test(itest).oddsRightPct mouse(imouse).expt(iexp).oddsRightPct];
            mouse(imouse).test(itest).testCon = [mouse(imouse).test(itest).testCon mouse(imouse).expt(iexp).pos(test_pos).intensities];
            mouse(imouse).test(itest).intensities = [mouse(imouse).test(itest).intensities; mouse(imouse).expt(iexp).intensities];
            mouse(imouse).posN = [mouse(imouse).posN; mouse(imouse).expt(iexp).posN];
            mouse(imouse).oddsRightPctN = [mouse(imouse).oddsRightPctN, mouse(imouse).expt(iexp).oddsRightPct];
            for ipos = 1:2
                mouse(imouse).test(itest).pos(ipos).nCorrect = [mouse(imouse).test(itest).pos(ipos).nCorrect; mouse(imouse).expt(iexp).pos(ipos).nCorrect];
                mouse(imouse).test(itest).pos(ipos).nIncorrect = [mouse(imouse).test(itest).pos(ipos).nIncorrect; mouse(imouse).expt(iexp).pos(ipos).nIncorrect];
            end
        end
        for itest = 1:3
            for ipos = 1:2
                if length(mouse(imouse).test(itest).pos(ipos).nCorrect + mouse(imouse).test(itest).pos(ipos).nIncorrect)>1
                    [x y] = binofit(sum(mouse(imouse).test(itest).pos(ipos).nCorrect,1),sum(mouse(imouse).test(itest).pos(ipos).nCorrect,1) + sum(mouse(imouse).test(itest).pos(ipos).nIncorrect,1));
                    mouse(imouse).test(itest).pos(ipos).pctCorrect = x;
                    mouse(imouse).test(itest).pos(ipos).ci95 = y;
                end
            end
        end
    end
    save([rc.fitOutputSummary '\' date '-2afc-attn_mouse.mat'], 'mouse')
end

    
    