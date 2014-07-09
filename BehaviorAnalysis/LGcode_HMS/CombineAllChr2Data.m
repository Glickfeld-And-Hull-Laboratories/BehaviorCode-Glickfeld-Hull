function CombineAllChr2Data(itask)

    mouse = createExptStructPV(itask);
    pv = behavParamsPV;
    rc = behavConstsHADC8;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

    for imouse = 1:length(pv.mouse_mat)
        if pv.chr2_mat(imouse) == 1
            ipos = 1;
            if ~isempty(mouse(imouse).task(itask).pos(ipos).ind)
                for ipow = 1:length(pv.power_mat)
                    pow_exist(ipow) = ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind);
                end
                pow_test = sum(pow_exist);
                start = 1;
                figure;
                for ipow = 1:length(pv.power_mat)
                    if ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                        subplot(ceil(sqrt(pow_test)), ceil(sqrt(pow_test)), start)
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).intensities = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).success = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).early = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).missed = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).b2ind = [];
                        mouse(imouse).task(itask).pos(ipos).pow(ipow).whichtrials = [];
                        days = 0;
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                            if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh)
                                if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh)
                                    days = days+1;
                                    ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp);
                                    fn = fullfile('/Users/lindsey/Desktop/Data', ['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(ind)) '.mat']);
                                    load(fn);
                                    if itask == 1
                                        xval = cell2mat(input.gratingContrast);
                                        doContrast = 1;
                                        xval_sub = bsxfun(@minus,xval,mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh);
                                        xval_sub(find(xval_sub>=.1)) = chop(xval_sub(find(xval_sub>=.1)), 1);
                                        xval_sub(find(xval_sub<=-.1)) = chop(xval_sub(find(xval_sub<=-.1)), 1);
                                        xval_sub(find(xval_sub<.1 & xval_sub>-.1)) = roundto(xval_sub(find(xval_sub<.1 & xval_sub>-.1)),1);
                                        xval_norm = xval_sub;
                                    elseif itask == 2
                                        xval = cell2mat(input.gratingDirectionDeg);
                                        doOri = 1;
                                        xval_sub = bsxfun(@minus,xval,mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh);
                                        xval_sub(find(xval_sub>=10)) = chop(xval_sub(find(xval_sub>=10)), 2, 4);
                                        xval_sub(find(xval_sub<=-10)) = chop(xval_sub(find(xval_sub<=-10)), 2, 4);
                                        xval_sub(find(xval_sub>-10 & xval_sub<10)) = chop(xval_sub(find(xval_sub>-10 & xval_sub<10)), 1, 4);   
                                        xval_sub(find(xval_sub<1 & xval_sub>-1)) = 0;
                                        xval_norm = xval_sub;
                                    end
                                    trials = cell2mat(xd.TrialRangeToUse(ind));
                                    if isempty(trials)
                                        trials = input.trialSinceReset;
                                        range = 1:trials;
                                    else
                                        range = str2num(trials);
                                    end
                                    trial_mat = zeros(size(xval));
                                    trial_mat(1,range) = 1;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).intensities = [mouse(imouse).task(itask).pos(ipos).pow(ipow).intensities xval_norm];
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).success = [mouse(imouse).task(itask).pos(ipos).pow(ipow).success strcmp(input.trialOutcomeCell, 'success')];
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).early = [mouse(imouse).task(itask).pos(ipos).pow(ipow).early strcmp(input.trialOutcomeCell, 'failure')];
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).missed = [mouse(imouse).task(itask).pos(ipos).pow(ipow).missed strcmp(input.trialOutcomeCell, 'ignore')];                    
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).b2ind = [mouse(imouse).task(itask).pos(ipos).pow(ipow).b2ind cell2mat_padded(input.tBlock2TrialNumber)];
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).whichtrials = [mouse(imouse).task(itask).pos(ipos).pow(ipow).whichtrials trial_mat];
                                end
                            end
                        end
                        block2V = mouse(imouse).task(itask).pos(ipos).pow(ipow).b2ind;
                        intensityV = mouse(imouse).task(itask).pos(ipos).pow(ipow).intensities;
                        successIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).success;
                        earlyIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).early;
                        missedIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).missed;
                        whichtrials = mouse(imouse).task(itask).pos(ipos).pow(ipow).whichtrials;
                        if all(isnan(block2V))
                            block2V = zeros(size(block2V));
                        end
                        block2Indices = unique(block2V(~isnan(block2V)));
                        nBlock2Indices = length(block2Indices);
                        assert(nBlock2Indices > 0 && nBlock2Indices <= 2);
                        for iB = 1:nBlock2Indices
                            tBN = block2Indices(iB);
                            b2Ix = block2V == tBN;

                            intensitiesC{iB} = unique(intensityV(b2Ix));
                            nIntensities(iB) = length(intensitiesC{iB});

                            intensityNums = block2V*NaN;
                            for iI = 1:nIntensities(iB)
                                tI = intensitiesC{iB}(iI);
                                iIx = intensityV == tI;

                                nCorr(iB, iI) = sum(iIx & successIx & b2Ix);
                                nEarly(iB, iI) = sum(iIx & earlyIx & b2Ix);
                                nMiss(iB, iI) = sum(iIx & missedIx & b2Ix);
                                nRawTot(iB, iI) = sum(iIx & b2Ix);
                            end
                            nCorrPlusMissC{iB} = nCorr(iB,1:nIntensities(iB)) + nMiss(iB,1:nIntensities(iB));
                        end
                        nCorrPlusMiss = nCorr+nMiss;
                        pctCorr = nCorr./nCorrPlusMiss;
                        pctCorr(pctCorr==1) = pctCorr(pctCorr==1)-10*eps;
                        percentsCorrect = pctCorr;
                        plot(intensitiesC{1},percentsCorrect(1,1:length(intensitiesC{1})),'ok')
                        hold on; 
                        plot(intensitiesC{2},percentsCorrect(2,1:length(intensitiesC{2})),'or')
                        if itask == 1
                            xlim([-.5 1]);
                        elseif itask == 2
                            xlim([-20 50]);
                        end
                        title(['Power = ' num2str(pv.power_mat(ipow)) '; n = ' num2str(days) 'days']);
                        start = start+1;
                    end
                end
                suptitle(['Summary for mouse i' num2str(pv.mouse_mat(imouse))])
                pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-combined-summary-i' num2str(pv.mouse_mat(imouse)) '-' date '.pdf']);
                exportfig_print(gcf, pn, 'FileFormat', 'pdf');
            end
        end
    end
end
%             fieldNames = { 'intensitiesC', 'percentsCorrect', ...
%         'nIntensities', 'nCorrPlusMissC', ...
%         'doContrast', 'doOri', ...
%         'successIx', 'earlyIx', 'missedIx', 'intensityV', ...
%         'block2V', 'block2Indices', 'nBlock2Indices', ...
%         'nMiss', 'nCorr', 'nCorrPlusMiss', 'whichTrials' };
% 
%             nFs = length(fieldNames);
%             for iF=1:nFs
%                 tFN = fieldNames{iF};
%                 bs.(tFN) = eval([tFN ';']);
%             end
% 
%         end
% end
