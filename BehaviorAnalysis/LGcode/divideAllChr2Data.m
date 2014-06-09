function divideAllChr2Data(itask)

    mouse = createExptStructPV(itask);
    pv = behavParamsPV;
    rc = behavConstsHADC8;
    xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);

    for imouse = 1
        if pv.chr2_mat(imouse) == 1
            ipos = 1;
            if ~isempty(mouse(imouse).task(itask).pos(ipos).ind)
                for ipow = 1:length(pv.power_mat)
                    pow_exist(ipow) = ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind);
                end
                pow_test = sum(pow_exist);
                for ipow = 1:length(pv.power_mat)
                    if ~isempty(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                        for iexp = 1:length(mouse(imouse).task(itask).pos(ipos).pow(ipow).ind)
                            if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(1).thresh)
                                if ~isnan(mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).laser(2).thresh)
                                    figure;
                                    ind = mouse(imouse).task(itask).pos(ipos).pow(ipow).ind(iexp);
                                    fn = fullfile('/Users/lindsey/Desktop/Data', ['data-i' num2str(pv.mouse_mat(imouse)) '-' cell2mat(xd.DateStr(ind)) '.mat']);
                                    load(fn);
                                    if itask == 1
                                        xval = cell2mat(input.gratingContrast);
                                        doContrast = 1;
                                    elseif itask == 2
                                        xval = cell2mat(input.gratingDirectionDeg);
                                        doOri = 1;
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
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).intensities = xval;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).success = strcmp(input.trialOutcomeCell, 'success');
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early = strcmp(input.trialOutcomeCell, 'failure');
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).missed = strcmp(input.trialOutcomeCell, 'ignore');                    
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).b2ind = cell2mat_padded(input.tBlock2TrialNumber);
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).whichtrials = trial_mat;
                                    mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).firsthalftrials = cell2mat_padded(input.reqHoldTimeMs)<ceil(max(cell2mat_padded(input.reqHoldTimeMs))./2);
                                end
                            end
                            block2V = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).b2ind;
                            intensityV = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).intensities;
                            successIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).success;
                            earlyIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).early;
                            missedIx = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).missed;
                            whichtrials = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).whichtrials;
                            firsthalftrials = mouse(imouse).task(itask).pos(ipos).pow(ipow).expt(iexp).firsthalftrials;
                            if all(isnan(block2V))
                                block2V = zeros(size(block2V));
                            end
                            block2Indices = unique(block2V(~isnan(block2V)));
                            nBlock2Indices = length(block2Indices);
                            assert(nBlock2Indices > 0 && nBlock2Indices <= 2);
                            halfIndices = unique(firsthalftrials(~isnan(firsthalftrials)));
                            nHalfIndices = length(halfIndices);
                            for iB = 1:nBlock2Indices
                                tBN = block2Indices(iB);
                                b2Ix = block2V == tBN;

                                intensitiesC{iB} = unique(intensityV(b2Ix));
                                nIntensities(iB) = length(intensitiesC{iB});

                                intensityNums = block2V*NaN;
                                for iH = 1:2
                                    tH = halfIndices(iH);
                                    h2Ix = firsthalftrials == tH;
                                    for iI = 1:nIntensities(iB)
                                        tI = intensitiesC{iB}(iI);
                                        iIx = intensityV == tI;

                                        nCorr(iB, iI, iH) = sum(iIx & successIx & b2Ix & h2Ix);
                                        nEarly(iB, iI, iH) = sum(iIx & earlyIx & b2Ix & h2Ix);
                                        nMiss(iB, iI, iH) = sum(iIx & missedIx & b2Ix& h2Ix);
                                        nRawTot(iB, iI, iH) = sum(iIx & b2Ix & h2Ix);
                                    end
                                nCorrPlusMissC{iB,iH} = nCorr(iB,1:nIntensities(iB),iH) + nMiss(iB,1:nIntensities(iB),iH);
                                end
                                nCorrPlusMiss = nCorr+nMiss;
                                pctCorr = nCorr./nCorrPlusMiss;
                                pctCorr(pctCorr==1) = pctCorr(pctCorr==1)-10*eps;
                                percentsCorrect = pctCorr;
                            end
                            subplot(2,1,1)
                            plot(intensitiesC{1},percentsCorrect(1,1:length(intensitiesC{1}),1),'ok')
                            hold on; 
                            plot(intensitiesC{2},percentsCorrect(2,1:length(intensitiesC{2}),1),'or')
                            subplot(2,1,2)
                            plot(intensitiesC{1},percentsCorrect(1,1:length(intensitiesC{1}),2),'ok')
                            hold on; 
                            plot(intensitiesC{2},percentsCorrect(2,1:length(intensitiesC{2}),2),'or')
                            suptitle(['Mouse = i' num2str(pv.mouse_mat(imouse)) ' Power = ' num2str(pv.power_mat(ipow)) ' ' cell2mat(xd.DateStr(ind))]);
                        end
                    end
                end
%                 pn = fullfile(rc.fitOutputSummary, ['Chr2-Task' num2str(itask) '-combined-summary-i' num2str(pv.mouse_mat(imouse)) '-' date '.pdf']);
%                 exportfig_print(gcf, pn, 'FileFormat', 'pdf');
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
