clear all
close all

%%
%contrast discrim 
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\BehaviorSummaries';
load(fullfile(data_path, 'CurrMouseList.mat')) 
behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\';
spSz = {4,2};
save_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\BehaviorSummaries\';
info = struct();
sumTrials = zeros();
yColor = [1 1 0]*0.5;

%%


for m = 1:length(allmice)
    
    clearvars -except data_path behav_path spSz save_path info m allmice yColor...
        doDiscrim doOri doMask doAdapt doFeedback doBlock2 Reward FBMotion  incTO ignTO StatPeriod Iti LV1  LV2 LV3 LV4...
        reactTime stimOnTime StimContrast StimDiameter StimEccentricity mouse_name task_name fractMaskTrials RandReward
    
%     if m == 1
%         continue
%     end
    
    mouse = allmice(m).name
    task = allmice(m).task
    info(m).mouse = mouse;
    %% dates
   
   dates = [220628; 220629; 220630];
   folder_dates = dates;
   trials = cell(1, length(dates));
%    mat_use{1} = 2;
%    mat_use{2} = 2; 
%    mat_use{3} = 2;
%    mat_use{4} = 2;
%    mat_use{5} = 2;

   mat_use{1} = [];
   mat_use{2} = []; 
   mat_use{3} = [];
   mat_use{4} = [];
   mat_use{5} = [];
   mat_use{6} = [];

   fprintf([mouse '\n']);
   
   dates_temp = dates;
   
   % check for mat files on dates
   for id = 1:size(dates,1)
        check_expt_mat = dir(fullfile(behav_path, ['data-' mouse '-' num2str(dates(id,:)) '-*']));
        if size(check_expt_mat,1) == 0
            rmidx = find(dates_temp == dates(id));
            dates_temp(rmidx) = [];        
        end
   end
   
   dates = dates_temp;

    %% variable for color scale
    cs = fliplr(0:.90/(length(dates)-1):.90);
    cmap = flipud(colormap(gray(length(dates)+1)));
    cmap(1,:) = [];
    
    cmap2 = flipud(colormap(winter(length(dates)+1)));
    cmap2(1,:) = [];
    
    %% 
    clear temp;
    figure;
    for id = 1:size(dates,1)
        info(m).dates = dates;
        expt_mat = dir(fullfile(behav_path, ['data-' mouse '-' num2str(dates(id,:)) '-*']));
        if size(expt_mat,1) > 1   
            if ~isempty(mat_use{id})
                load(fullfile(behav_path,expt_mat(mat_use{id}).name))
            else 
%                 x = zeros(1,size(expt_mat,1));
                matsize = [];
                for im = 1:size(expt_mat,1)
                    load(fullfile(behav_path,expt_mat(im).name))
                    if isfield(input,'trialOutcomeCell')
                        matsize(im) = expt_mat(im).bytes;           
                    end
                end
                biggest_mat_id = find(matsize == max(matsize));
                load(fullfile(behav_path,expt_mat(biggest_mat_id).name))
%                         x(1,im) = 1;
%                     end
                
%                 if sum(x,2) == 1
%                     load(fullfile(behav_path,expt_mat(find(x==1)).name))
% 
%                 elseif sum(x,2)>1
%                     error('Too many mat files')
%                 elseif sum(x,2)==0
%                     error('No mat files')
%                 end
            end
        else
            load(fullfile(behav_path,expt_mat.name))
        end
%             [s, b] = selectCalc(input,trials{id});
            [s, ms] = selectCalcLikeLego(input);
            select(id) = s;
            info(m).select = select;
            
            if input.doMask
%                 [s2, b2] = selectCalc2(input,trials{id});
                select2(id) = ms;
                info(m).select2 = select2;
            end
            
                if ~isempty(trials{id})
                    input = trialChopper(input,trials{id});
                end
                
                block2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
                block1Ix = block2Ix==0;
                correctIx = strcmp(input.trialOutcomeCell, 'success');
                incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
                ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
                missIx = ignoreIx;
                if isfield(input,'doMask')
                    maskIx = celleqel2mat_padded(input.tDoMask);
                    type1MaskIx = maskIx;
                    if isfield(input,'doType2Mask')
                      type2MaskIx = celleqel2mat_padded(input.tDoType2Mask);
                      type1MaskIx(find(type2MaskIx)) = 0;
                    else
                      type2MaskIx = zeros(size(block2Ix));
                    end
                else
                    maskIx = zeros(size(block2Ix));
                    type2MaskIx = zeros(size(block2Ix));
                    type1MaskIx = zeros(size(block2Ix));
                end
                nCorr(id) = sum(correctIx);
                nInc(id) = sum(incorrectIx);
                nIg(id) = sum(ignoreIx);
                nTrials(id) = length(input.trialOutcomeCell);

                leftTrialIx = celleqel2mat_padded(input.tLeftTrial);
                leftTrNs = find(leftTrialIx);
                rightTrialIx = ~leftTrialIx;
                rightTrNs = find(rightTrialIx);
                leftTrialIx = logical(leftTrialIx);
                nLeft(id) = sum(leftTrialIx);
                nRight(id) = sum(rightTrialIx);

                leftOutcomes = input.trialOutcomeCell(leftTrialIx);
                leftCorr = strcmp(leftOutcomes, 'success');
                nLeftCorr(id) = sum(leftCorr);
                leftIgn = strcmp(leftOutcomes, 'ignore');
                nLeftIgn(id) = sum(leftIgn);
                leftInc = strcmp(leftOutcomes, 'incorrect');
                nLeftInc(id) = sum(leftInc);

                rightOutcomes = input.trialOutcomeCell(rightTrialIx);
                rightCorr = strcmp(rightOutcomes, 'success');
                nRightCorr(id) = sum(rightCorr);
                rightIgn = strcmp(rightOutcomes, 'ignore');
                nRightIgn(id) = sum(rightIgn);
                rightInc = strcmp(rightOutcomes, 'incorrect');
                nRightInc(id) = sum(rightInc);  
                
                
                %% decision time cdf
                
                axH = subplot(spSz{:}, 2);
                hold on
                decisionMax = ceil(input.reactionTimeMs/1000)*1000;
                if decisionMax <= 6000
                    decisionInterval = 1000;
                else
                    decisionInterval = 2000;
                end
                if decisionMax <= 10000
                    decMax = decisionMax;
                else
                    decMax = 10000;
                end

                if find(isnan(block2Ix))
                  block2Ix(find(isnan(block2Ix))) = 0;
                end

                if sum(block2Ix)== 0 % block 2 off
                  if nLeft>0
                    pH1 = cdfplot([input.tDecisionTimeMs{leftTrialIx}]);
                      set(pH1, 'Color', 'b');
                  end
                  if nRight>0
                    pH2 = cdfplot([input.tDecisionTimeMs{rightTrialIx}]);
                      set(pH2, 'Color', [0.8 0.8 0]);
                  set(gca, 'XLim', [0 decMax], ...
                          'YLim', [0 1], ...
                          'XTick', [0:decisionInterval:decMax]);
                  end
                elseif sum(block2Ix)>0 % block 2 on
                  pH1 = cdfplot([input.tDecisionTimeMs{block2Ix == 0 & leftTrialIx == 1}]);
                    set(pH1,'Color', 'b', 'LineWidth', 2);
                  hold on
                  pH2 = cdfplot([input.tDecisionTimeMs{block2Ix == 0 & leftTrialIx == 0}]);
                    set(pH2, 'Color', [0.8 0.8 0],  'LineWidth', 2);
                  pH3 = cdfplot([input.tDecisionTimeMs{block2Ix == 1 & leftTrialIx == 1}]);
                    set(pH3, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
                  pH4 = cdfplot([input.tDecisionTimeMs{block2Ix == 1 & leftTrialIx == 0}]);
                    set(pH4, 'Color', [0.8 0.8 0], 'LineStyle', '--', 'LineWidth', 2);
                    set(gca, 'XLim', [0 decMax], ...
                          'YLim', [0 1], ...
                          'XTick', [0:decisionInterval:decMax]);
                end
                grid on;
                hold on;
                title('Decision Time CDF: y=right,b=left');
                xlabel('Time');
                
                %% Right/Left Contrast Decision
                
               axH = subplot(spSz{:}, 5);
               hold on
                if isfield(input,'doContrastDiscrim')
                  if input.doContrastDiscrim
                    if input.doContrastDetect
                      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
                    elseif input.gratingContrastDiffSPO > 10
                      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
                    elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
                      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
                    elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
                      differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
                    end
                  elseif input.doOriDiscrim
                    differenceRight = chop(celleqel2mat_padded(input.tGratingDirectionStart) - double(input.gratingTargetDirection),2);
                    if isfield(input,'doMask')
                      differenceRight_mask = chop(celleqel2mat_padded(input.tPlaidDirectionStart) - double(input.gratingTargetDirection),2);
                    else
                      differenceRight_mask = nan(size(differenceRight));
                    end
                  end
                else
                  if input.gratingContrastDiffSPO > 10
                    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
                  elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
                    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
                  elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
                    differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
                  end
                end
                

                plotTrsB1 = differenceRight((correctIx|incorrectIx)&~block2Ix&~maskIx);
                nLevelsB1 = unique(plotTrsB1);
                percentContCellB1 = cell(1,length(nLevelsB1));
                for kk=1:length(nLevelsB1)
                    valB1 = nLevelsB1(kk);
                    valIxB1 = differenceRight==valB1;
                    totalNTrialsValB1 = length(differenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix&~maskIx));
                    if min(differenceRight) < 0
                      if valB1>=0,
                          ind = setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]);  
                          rightNTrialsValB1 = length(ind);
                          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
                      elseif valB1<0,
                          ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]);
                          rightNTrialsValB1 = length(ind);
                          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
                      end
                    else
                        if valB1>=1,
                            ind = setdiff(intersect(find(valIxB1),find(correctIx)), find(block2Ix));
                            rightNTrialsValB1 = length(ind);
                            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
                        elseif valB1<1,
                            ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), find(block2Ix));
                            rightNTrialsValB1 = length(ind);
                            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
                        end
                    end
                end

                if sum(maskIx)>0
                    if sum(type1MaskIx)>0
                      plotTrsM1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type1MaskIx);
                      nLevelsM1 = unique(plotTrsM1);
                      percentContCellM1 = cell(1,length(nLevelsM1));
                      for kk=1:length(nLevelsM1)
                          valM1 = nLevelsM1(kk);
                          valIxM1 = differenceRight_mask==valM1;
                          totalNTrialsValM1 = length(differenceRight_mask(valIxM1&(correctIx|incorrectIx)&~block2Ix&type1MaskIx));
                          if min(differenceRight_mask) < 0
                            if valM1>=0,
                                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                                rightNTrialsValM1 = length(ind);
                                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
                            elseif valM1<0,
                                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                                rightNTrialsValM1 = length(ind);
                                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
                            end
                          else
                              if valM1>=1,
                                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                                  rightNTrialsValM1 = length(ind);
                                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
                              elseif valM1<1,
                                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                                  rightNTrialsValM1 = length(ind);
                                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
                              end
                          end
                      end
                    end
                    if sum(type2MaskIx)>0
                      plotTrsT2M1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type2MaskIx);
                      nLevelsT2M1 = unique(plotTrsT2M1);
                      percentContCellT2M1 = cell(1,length(nLevelsT2M1));
                      for kk=1:length(nLevelsT2M1)
                          valT2M1 = nLevelsT2M1(kk);
                          valIxT2M1 = differenceRight_mask==valT2M1;
                          totalNTrialsValT2M1 = length(differenceRight_mask(valIxT2M1&(correctIx|incorrectIx)&~block2Ix&type2MaskIx));
                          if min(differenceRight_mask) < 0
                            if valT2M1>=0,
                                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                                rightNTrialsValT2M1 = length(ind);
                                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
                            elseif valT2M1<0,
                                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                                rightNTrialsValT2M1 = length(ind);
                                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
                            end
                          else
                              if valT2M1>=1,
                                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                                  rightNTrialsValT2M1 = length(ind);
                                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
                              elseif valT2M1<1,
                                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                                  rightNTrialsValT2M1 = length(ind);
                                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
                              end
                          end
                      end
                    end
                end

                if sum(block2Ix)>0
                  plotTrsB2 = differenceRight((correctIx|incorrectIx)&block2Ix);
                  nLevelsB2 = unique(plotTrsB2);
                  percentContCellB2 = cell(1,length(nLevelsB2));
                  for kk=1:length(nLevelsB2)
                    valB2 = nLevelsB2(kk);
                    valIxB2 = differenceRight==valB2;
                    totalNTrialsValB2 = length(differenceRight(valIxB2&(correctIx|incorrectIx)&block2Ix));
                    if min(differenceRight) < 0
                      if valB2>=0,
                          ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
                          rightNTrialsValB2 = length(ind);
                          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
                      elseif valB2<0,
                          ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
                          rightNTrialsValB2 = length(ind);
                          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
                      end
                    else
                        if valB2>=1,
                            ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
                            rightNTrialsValB2 = length(ind);
                            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
                        elseif valB2<1,
                            ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
                            rightNTrialsValB2 = length(ind);
                            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
                        end
                    end
                  end
                end

                if sum(maskIx&block2Ix)>0
                    if sum(type1MaskIx)>0 type1MaskIx
                        plotTrsM2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type1MaskIx);
                        nLevelsM2 = unique(plotTrsM2);
                        percentContCellM2 = cell(1,length(nLevelsM2));
                        for kk=1:length(nLevelsM2)
                            valM2 = nLevelsM2(kk);
                            valIxM2 = differenceRight_mask==valM2;
                            totalNTrialsValM2 = length(differenceRight_mask(valIxM2&(correctIx|incorrectIx)&block2Ix&type1MaskIx));
                            if min(differenceRight_mask) < 0
                              if valM2>=0,
                                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
                                  rightNTrialsValM2 = length(ind);
                                  percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
                              elseif valM2<0,
                                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
                                  rightNTrialsValM2 = length(ind);
                                  percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
                              end
                            else
                                if valM2>=1,
                                    ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
                                    rightNTrialsValM2 = length(ind);
                                    percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
                                elseif valM2<1,
                                    ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
                                    rightNTrialsValM2 = length(ind);
                                    percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
                                end
                            end
                        end
                     end
                      if sum(type2MaskIx)>0
                        plotTrsT2M2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type2MaskIx);
                        nLevelsT2M2 = unique(plotTrsT2M2);
                        percentContCellT2M2 = cell(1,length(nLevelsT2M2));
                        for kk=1:length(nLevelsT2M2)
                            valT2M2 = nLevelsT2M2(kk);
                            valIxT2M2 = differenceRight_mask==valT2M2;
                            totalNTrialsValT2M2 = length(differenceRight_mask(valIxT2M2&(correctIx|incorrectIx)&block2Ix&type2MaskIx));
                            if min(differenceRight_mask) < 0
                              if valT2M2>=0,
                                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
                                  rightNTrialsValT2M2 = length(ind);
                                  percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
                              elseif valT2M2<0,
                                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
                                  rightNTrialsValT2M2 = length(ind);
                                  percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
                              end
                            else
                                if valT2M2>=1,
                                    ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
                                    rightNTrialsValT2M2 = length(ind);
                                    percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
                                elseif valT2M2<1,
                                    ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
                                    rightNTrialsValT2M2 = length(ind);
                                    percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
                                end
                            end
                        end
                      end
                 end
                

                if min(differenceRight) < 0
                    minX = min(differenceRight);
                    maxX = max(differenceRight);
                    xLimm = [minX maxX];
                else
                    minX = min(differenceRight,[],2);
                    maxX = max(differenceRight,[],2);
                    xLimm = [minX maxX];
                        if ~(xLimm(1)==0),
                            xL1 = [floor(log10(xLimm(1))) ceil(log10(xLimm(2)))];
                        else
                            xL1 = [0 ceil(log10(xLimm(2)))];
                        end
                    xTickL = 10.^(xL1(1):1:xL1(2));
                    xTickL = xTickL(xTickL>=xLimm(1) & xTickL<=xLimm(2));
                    xTLabelL = cellstr(num2str(xTickL(:)));
                end

                nLevelsB1 = nLevelsB1(~isnan(nLevelsB1));
                pH1 = plot(nLevelsB1, cell2mat(percentContCellB1), 'Color', [cs(id) cs(id) cs(id)], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                if sum(block2Ix)>= 1
                  pH2 = plot(nLevelsB2, cell2mat(percentContCellB2), 'Color', [cs(id) cs(id) 1], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                end

                if sum(type1MaskIx)>= 1
                  pH4 = plot(nLevelsM1, cell2mat(percentContCellM1), 'Color', [cs(id) cs(id) 1], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                end
                if sum(type1MaskIx&block2Ix)>= 1
                  pH5 = plot(nLevelsM2, cell2mat(percentContCellM2),  'Linestyle', '--', 'Color', [cs(id) cs(id) 1], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                end
                
                if sum(type2MaskIx)>= 1
                  pH6 = plot(nLevelsT2M1, cell2mat(percentContCellT2M1), 'Color', [cs(id) 1 cs(id)], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                end
                if sum(type2MaskIx&block2Ix)>= 1
                  pH7 = plot(nLevelsT2M2, cell2mat(percentContCellT2M2), 'Linestyle', '--', 'Color', [cs(id) 1 cs(id)], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 8);
                end


                if input.gratingContrastDiffSPO <= 100
                  vH = plot([1 1],[0 1]);
                else
                  vH = plot([0 0],[0 1]);
                end
                set(vH, 'Color', 'g');

                set(gca, 'XLim', [minX maxX], ...
                         'YLim', [0 1]);
                if min(differenceRight) < 0
                  if isfield(input,'doContrastDiscrim')
                    if input.doContrastDiscrim
                        xlabel('Contrast Difference (R-L)')
                        set(gca, 'XTick', [-1:0.5:1], ...
                                 'YTick', [0:0.25:1],...
                                 'XGrid', 'on');
                    elseif input.doOriDiscrim
                        xlabel('Ori Difference')
                        set(gca, 'XTick', [-input.gratingMaxDirectionDiff:15:input.gratingMaxDirectionDiff], ...
                                 'YTick', [0:0.25:1],...
                                 'XGrid', 'on');
                    end
                  else
                    xlabel('Contrast Difference (R-L)')
                        set(gca, 'XTick', [-1:0.5:1], ...
                                 'YTick', [0:0.25:1],...
                                 'XGrid', 'on');
                  end
                else
                  set(gca,'XScale', 'log', ...
                            'XGrid', 'on',...
                            'XTick', xTickL,...
                            'XTickLabel', xTLabelL);
                  if isfield(input,'doContrastDiscrim')
                    if input.doContrastDiscrim
                      xlabel('Contrast Difference (R/L)')
                    elseif input.doOriDiscrim
                      xlabel('Ori Difference')
                    end
                  else
                    xlabel('Contrast Difference (R/L)')
                  end
                end

                if input.doOriDiscrim
                    ylabel('% Left')
                else
                    ylabel('% Right')
                end
                grid on

  rawDT = celleqel2mat_padded(input.tDecisionTimeMs);
    MeanDecisionTime(id) =  mean(rawDT(correctIx|incorrectIx));
    info(m).meanDT = MeanDecisionTime;
                
               
    end
    

    
    
%     title({['s = ' num2str(chop(select,2))]})
      title('% choice by stimulus')
        % nTrials left and right

%             days = [1 2 3];

    % nLeft trials by day
    
    x = 1:length(dates);


    axH(1) = subplot(spSz{:}, 3);
    plot(x, nLeftCorr, 'Color', 'k', 'LineWidth', 2)
    hold on
    plot(x, nLeftInc, 'Color', 'g', 'LineWidth', 2)
    plot(x, nLeftIgn, 'Color', 'm', 'LineWidth', 2)
    plot(x, nLeft,  '--', 'Color', 'k', 'LineWidth', 2)
    title('Left')
    ylabel('nTrials')
    xlabel('Days')

    xticks(x)
    for g = 1:length(x)
        datestring = num2str(dates(g));
        day{g} = datestring(3:6);
    end
    xticklabels(day)

    % nRight trials by day

    axH(2) = subplot(spSz{:}, 4);
    plot(x, nRightCorr, 'Color', 'k', 'LineWidth', 2)
    hold on
    plot(x, nRightInc, 'Color', 'g', 'LineWidth', 2)
    plot(x, nRightIgn, 'Color', 'm', 'LineWidth', 2)
    plot(x, nRight, '--', 'Color', 'k', 'LineWidth', 2)
    title('Right')
    ylabel('nTrials')
    xlabel('Days')

    xticks(x)
    for g = 1:length(x)
        datestring = num2str(dates(g));
        day{g} = datestring(3:6);
    end
    xticklabels(day)
    
    linkaxes(axH, 'y');

   % n day performance
    
    axH(3) = subplot(spSz{:}, 6);
    plot(x, select, 'Color', 'k', 'LineWidth', 2)
    hold on
    if input.doMask
        plot(x, select2, 'Color', 'b', 'LineWidth', 2)
    end
    hline(0.9)
    hline(0)
    title([ num2str(length(x)) ' Day performance'])
    ylabel('Selectivity')
    xlabel('Days')
    xticks(x)
    ylim([-0.1 1])
    for g = 1:length(x)
        datestring = num2str(dates(g));
        day{g} = datestring(3:6);
    end
    xticklabels(day)
    
    yyaxis right
    plot(x, MeanDecisionTime, 'LineWidth', 2)
    ylabel('Decision Time')
    
   

    %% performance values - totals and current settings (make work)
    
    
        total_nCorr = sum(nCorr);
        total_nInc = sum(nInc);
        total_nIg = sum(nIg);
        total_nTrials = sum(nTrials);
        total_nLeft = sum(nLeft);
        total_nRight = sum(nRight);
        total_nLeftCorr = sum(nLeftCorr);
        total_nLeftIgn = sum(nLeftIgn);
        total_nLeftInc = sum(nLeftInc);
        total_nRightCorr = sum(nRightCorr);
        total_nRightIgn = sum(nRightIgn);
        total_nRightInc = sum(nRightInc);
        
        ndays = length(x);
        
        meanCorr = roundto(total_nCorr/ndays,0);
        meanInc = roundto(total_nInc/ndays,0);
        meanIg = roundto(total_nIg/ndays,0);
        meanTrials = roundto(total_nTrials/ndays,0);
        
        
        
        
        
        axI = subplot(spSz{:}, 1);						% default axes are 0 to 1
        set(axI, 'Visible', 'off');
        set(axI, 'OuterPosition', [0.02 0.75, 0.25, 0.2])       
%         set(gcf, 'Visible', 'off'); % hide figure during text
%                                     % drawing - kludge
        text(0.00, 1.25, [mouse ': ' task ' : ' [day{1} ' - ' day{end}]], 'FontWeight', 'bold', 'FontSize', 16);
        hold on
        t2H(1) = text(0.00, 1.15, {'    ', 'Trials:', 'Correct:', 'Incorrect:', 'Missed:'});
        t2H(2) = text(0.4, 1.15, {'Left', sprintf('%d', total_nLeft), sprintf('%d', total_nLeftCorr), ...
				sprintf('%d', total_nLeftInc), sprintf('%d', total_nLeftIgn)});
        t2H(3) = text(0.65, 1.15, {'Right', sprintf('%d', total_nRight), sprintf('%d', total_nRightCorr), ...
				sprintf('%d', total_nRightInc), sprintf('%d', total_nRightIgn)});
        t2H(4) = text(0.9, 1.15, {'Total', sprintf('%d', total_nTrials), sprintf('%d', total_nCorr), ...
                sprintf('%d', total_nInc), sprintf('%d', total_nIg)});
        t2H(5) = text(1.15, 1.15, {'avg/day', sprintf('%d', meanTrials), sprintf('%d', meanCorr), ...
                sprintf('%d', meanInc), sprintf('%d', meanIg)});
        set(t2H, 'VerticalAlignment', 'top', ...
                 'HorizontalAlignment', 'left');
             
             
             
        trPer80Str = regexprep(['[' num2str(input.trPer80V, '%3d') ']'], '[ *', '[');
        rewStr = mat2str(input.rewardTimeUs/1000);

             
        if input.doFeedbackMotion
            FB_str = 'On';
        else
            FB_str = 'Off';
        end
        
        if ~isfield(input, 'doRandReward')
            input.doRandReward = 0;
        end
        
        if input.doRandReward
            randRW_str = ' -rand';
        else
            randRW_str = '-norm';
        end
        
        size = num2str(min(input.gratingMaxDiameterDeg));
        ecc =  num2str(min(input.gratingEccentricityDeg));
             
         tStr = sprintf( ['Decision Time: \t%5.2f s;   ITI %d ms \n', ...
             'Stationary Period: \t%5d ms \n', ...
             ['Stim On Time: \t%5d ms; sz: ' size '; ecc: ' ecc ' \n'], ...
             ['Too-Fast Time: %d ms;  Reward: %s ms' randRW_str '\n'], ...
             'Timeouts (ign,inc):\t%4.1f, %4.1f s\n', ...
             'trPer80: %s\n', ...
             ['Feedback' FB_str ', Motion Sensitivity: \t%2.3f' '; LEDmw: \t%2.3f \n']], ...
            input.reactionTimeMs/1000, ...
            input.itiTimeMs, ...
            input.stationaryPeriodMs,...
            input.stimOnTimeMs, ...
            input.tooFastTimeMs, ...
            rewStr, ...
            input.ignoreTimeoutMs/1000, ...
            input.incorrectTimeoutMs/1000, ...
            trPer80Str, ...
            input.feedbackMotionSensitivity,...
            input.block2TrialLaserPowerMw);

        text(0.0, 0.55, tStr, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');
% 


    %% DT vs select 
    axI = subplot(spSz{:}, 7);	
    DT = info(m).meanDT;
    S = info(m).select;
    hold on
    sp = scatter(DT, S, 'filled')
    sp.CData = cmap;
    xlabel('Decision Time')
    ylabel('Selectivity')
    title('DT vs select')
    ylim([-0.1 1])
    vline(2000)
    hline(0.9)
    hline(0)
    
    
    if input.doMask
        S2 = info(m).select2;
        sp2 = scatter(DT, S2, 'filled')
        sp2.CData = cmap2;
    end

    %% Variables for variable summary table
    

    if ~isfield(input, 'doRandReward')
        input.doRandReward = 0;
    end

    mouse_name{m,:} = mouse;
    task_name{m,:} = task; 
    doDiscrim(m,:) = input.doContrastDiscrim;
    doOri(m,:) = input.doOriDiscrim;
    doMask(m,:) = input.doMask;
    doType2Mask(m,:) = input.doType2Mask;
    fractMaskTrials(m,:) = input.fractionMaskTrials;
    doBlock2(m,:) = input.doBlock2;
    doAdapt(m,:) = input.doAdapt;
    doFeedback(m,:) = input.doFeedbackMotion;
    Reward(m,:) = input.rewardTimeUs;
    RandReward(m,:) = input.doRandReward;
    FBMotion(m,:) = input.feedbackMotionSensitivity;
    incTO(m,:) = input.incorrectTimeoutMs;
    ignTO(m,:) = input.ignoreTimeoutMs;
    StatPeriod(m,:) = input.stationaryPeriodMs;
    Iti(m,:) = input.itiTimeMs;
    LV1(m,:) = input.trPer80Level1;
    LV2(m,:) = input.trPer80Level2;
    LV3(m,:) = input.trPer80Level3;
    LV4(m,:) = input.trPer80Level4;
    reactTime(m,:) = input.reactionTimeMs;
    stimOnTime(m,:) = input.stimOnTimeMs;
    StimContrast(m,:) = input.gratingMaxContrast;
    StimDiameter(m,:) = input.gratingMaxDiameterDeg;
    StimEccentricity(m,:)  = input.gratingEccentricityDeg;
    
    %% Training progress over time
    
%     info2 = struct();
%     expt_mat2 = dir(fullfile(behav_path, ['data-' mouse '-*']));
% %     if m == 9 % i476 swtiched tasks on jun 8 2020 -Lego notes
% %         expt_mat2 = expt_mat2(9:end);
% %     end
%     
% %     if length(expt_mat2) > 120 % limit to 3 months
% %         expt_mat2 = expt_mat2(end-119:end);
% %     end
% 
%     trials2 = cell(1, length(expt_mat2));
% %     tdays = 1:length(expt_mat2);
% 
%     for id2 = 1:length(expt_mat2)
%         if expt_mat2(id2).bytes == 0
%             continue
%         end
%         
%  
%         progress = [num2str(id2/length(expt_mat2) * 100) '%']
%         load(fullfile(behav_path,expt_mat2(id2).name))
%         if ~isfield(input, 'trialOutcomeCell') && ~isfield(input, 'dGratingContrast')
%             continue
%         end
%         [s, b] = selectCalc(input,trials2{id2});
% %         [s, ms] = selectCalcLikeLego(input);
%         select(id2) = s;
%     %     info2(m).select = select;
%     
% %         if input.doMask
% %             [s2, b2] = selectCalc2(input,trials2{id2});
% %             select2(id2) = s2;
% %         end
% %         
%     
%             if ~isempty(trials2{id2})
%                 input = trialChopper(input,trials2{id2});
%             end
%         correctIx = strcmp(input.trialOutcomeCell, 'success');
%         incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
%         rawDT = celleqel2mat_padded(input.tDecisionTimeMs);
%         MeanDecisionTime(id2) =  mean(rawDT(correctIx|incorrectIx));
%     %     info(m).meanDT = MeanDecisionTime;
% 
% 
%     end
% 
%     MeanDecisionTime(select == 0) = [];
%     select(select == 0) = [];
%     tdays = 1:length(select);
%     
    
%     [smooth_select, smooth_MeanDT, tdays] = getTrainingHistory(behav_path, mouse);
    [smooth_select, smooth_MeanDT, tdays, smooth_select2, tdays2] = getTrainingHistory2(behav_path, mouse);

    axI = subplot(spSz{:}, 8);	
    plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
    hold on
    ylim([-0.1 1])
    hline(0.9)
    hline(0)
    xlabel('Training Days')
    ylabel('Selectivity')
    title('Total Performance History (MA)')
    
    plot(tdays2, smooth_select2, 'b', 'LineWidth', 1.25);
    
    yyaxis right
    plot(tdays, smooth_MeanDT, 'LineWidth', 1.24)
    ylabel('Decision Time')
    
    
    
    
    %% final things



    

    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
    
    for g = 1:length(folder_dates)
    save_datestring = num2str(folder_dates(g));
    save_day{g} = save_datestring(3:6);
    end
    
%    uncomment below to save 
   	save_daterange = [save_day{1} '-' save_day{end}];
    if m ==1  
        mkdir(save_path, save_daterange)
    end
    cd(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\BehaviorSummaries\', save_daterange, '\'])
    print(mouse, '-dpdf')
    
 
        
 

end

% uncomment to save variable spreadsheet
T = table(mouse_name, task_name, doMask, doType2Mask, fractMaskTrials, doBlock2, doFeedback, stimOnTime, FBMotion, Reward, RandReward, incTO, ignTO, StatPeriod, Iti, LV1, LV2, LV3, LV4, reactTime, StimContrast, StimDiameter, StimEccentricity)
filename = 'variableinfo.xlsx';
writetable(T, filename)

%%  DT vs select (separate summary figure)

% spSz2 = {2,2};
% 
% axJ = subplot(spSz2{:}, 1);
% for m = 1:length(allmice)
%     
%     if m == 5
%         continue
%     end
%     
%     if strcmp(allmice(m).task, 'Plaid')
%         color = 'r';
%     elseif strcmp(allmice(m).task, 'OriDiscrim')
%         color = [0.4660 0.6740 0.1880];
%     elseif strcmp(allmice(m).task, 'ContrastDiscrim')
%         color = 'b';
%     end
%     
%     DT = info(m).meanDT;
%     S = info(m).select;
%     hold on
%     scatter(DT, S, 'filled', 'MarkerEdgeColor', color, 'MarkerFaceColor', color)
%     
%     
%     
% end
% 
% xlabel('Decision Time')
% ylabel('Selectivity')
% title('r=plaid,b=contrast,g=ori')
% ylim([-0.1 1])
% vline(2000)
% hline(0.9)
% hline(0)
% % lables = [{'Plaid'}, {'Contrast'}, {'Ori'}];
% % legend(labels, 'Location', 'southeastoutside')
% % 
% % clear labels
% 
% 
% % cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\')
% % print('all mice', '-dpdf')
% 
% 
% % % example mouse
% % figure;
% % DT = info(1).meanDT;
% % S = info(1).select;
% % hold on
% % scatter(DT, S)
% % xlabel('Decision Time')
% % ylabel('Selectivity')
% % title('i484, 5 days')
% % ylim([0 1])
% % hline(0.9)
% % vline(2000)
% 
% % cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\')
% % print('example_i484', '-dpdf')
% 
% % by experiment
% 
% axJ = subplot(spSz2{:}, 2);
% for m = 1:4 % plaid
%     
%     if m == 5
%         continue
%     end
%     
%     DT = info(m).meanDT;
%     S = info(m).select;
%     hold on
%     scatter(DT, S, 'filled')
%     
%     
%     labels{m} = allmice(m).name;
%        
% end
% 
% xlabel('Decision Time')
% ylabel('Selectivity')
% title('Plaid')
% ylim([-0.1 1])
% vline(2000)
% hline(0.9)
% hline(0)
% labels= labels(~cellfun('isempty',labels));
% 
% legend(labels, 'Location', 'southeastoutside')
% 
% clear labels
% 
% % cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\')
% % print('Plaid', '-dpdf')
% 
% axJ = subplot(spSz2{:}, 3);
% for m = [6 8 9 10 12]  % Ori
%     
%     if m == 5
%         continue
%     end
%     
%     DT = info(m).meanDT;
%     S = info(m).select;
%     hold on
%     scatter(DT, S, 'filled')
%     
%     
% 
%     labels{m} = allmice(m).name;
%        
% end
% 
% xlabel('Decision Time')
% ylabel('Selectivity')
% title('Ori Discrim')
% ylim([-0.1 1])
% vline(2000)
% hline(0.9)
% hline(0)
% labels= labels(~cellfun('isempty',labels));
% 
% legend(labels, 'Location', 'southeastoutside')
% 
% clear labels
% 
% 
% % cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\')
% % print('Ori', '-dpdf')
% 
% 
% axJ = subplot(spSz2{:}, 4);
% for m = [7 11 13]  % contrast
%     
%     if m == 5
%         continue
%     end
%     
%     DT = info(m).meanDT;
%     S = info(m).select;
%     hold on
%     scatter(DT, S, 'filled')
%     
% 
%     labels{m} = allmice(m).name;
%     
%        
% end
% 
% xlabel('Decision Time')
% ylabel('Selectivity')
% title('Contrast Discrim')
% ylim([-0.1 1])
% vline(2000)
% hline(0.9)
% hline(0)
% labels= labels(~cellfun('isempty',labels));
% legend(labels, 'Location', 'southeastoutside')
% 
% % cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\')
% % print('Contrast', '-dpdf')
% 
% % finalize figure
% 
% sgtitle(save_daterange)
% 
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])          

%uncomment to save
% cd(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\BehaviorSummaries\', save_daterange, '\'])
% print('Selectivity_v_DecisionTime', '-dpdf')





