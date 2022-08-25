function [smooth_select, smooth_MeanDT, tdays, smooth_select2, tdays2, stim_info] = getTrainingHistory2(behav_path, mouse)
%    behav_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\';
%     mouse = 'i484';
%     
    stim_info = struct();
    expt_mat2 = dir(fullfile(behav_path, ['data-' mouse '-*']));
%     if m == 9 % i476 swtiched tasks on jun 8 2020 -Lego notes
%         expt_mat2 = expt_mat2(9:end);
%     end
    
%     if length(expt_mat2) > 120 % limit to 3 months
%         expt_mat2 = expt_mat2(end-119:end);
%     end

    trials2 = cell(1, length(expt_mat2));
%     tdays = 1:length(expt_mat2);

    select = zeros(1,1000);
    select2 = zeros(1,1000);

    for id2 = 1:length(expt_mat2)
        
     
        
        if expt_mat2(id2).bytes <100000
            continue
        end
        
 
        progress = [num2str(id2/length(expt_mat2) * 100) '%']
        load(fullfile(behav_path,expt_mat2(id2).name))
        if ~isfield(input, 'trialOutcomeCell') && ~isfield(input, 'dGratingContrast')
            continue
        end
%         [s, b] = selectCalc(input,trials2{id2});
        [s, ms] = selectCalcLikeLego(input);
        
        
        
        
        select(id2) = s;
    %     info2(m).select = select;
    
        if input.doMask
%             [s2, b2] = selectCalc2(input,trials2{id2});
            select2(id2) = ms;
        end
%         
    
            if ~isempty(trials2{id2})
                input = trialChopper(input,trials2{id2});
            end
        correctIx = strcmp(input.trialOutcomeCell, 'success');
        incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
        rawDT = celleqel2mat_padded(input.tDecisionTimeMs);
        MeanDecisionTime(id2) =  mean(rawDT(correctIx|incorrectIx));
    %     info(m).meanDT = MeanDecisionTime;

        
        size(id2) = input.gratingMaxDiameterDeg;
        eccentricity(id2) = input.gratingEccentricityDeg;
        FB(id2) = input.doFeedbackMotion;
        FB_motion(id2) = input.feedbackMotionSensitivity;
        LVL1(id2) = input.trPer80Level1;
        LVL2(id2) = input.trPer80Level2;
        LVL3(id2) = input.trPer80Level3;
        LVL4(id2) = input.trPer80Level4;

    end

    select = select(1:find(select,1,'last'));
    select2 = select2(1:find(select,1,'last'));
    MeanDecisionTime(select == 0) = [];
    
    size(select == 0) = [];
    eccentricity(select == 0) = [];
    FB(select == 0) = [];
    FB_motion(select == 0) = [];
    LVL1(select == 0) = [];
    LVL2(select == 0) = [];
    LVL3(select == 0) = [];
    LVL4(select == 0) = [];
    
    stim_info.size = size;
    stim_info.eccentricity = eccentricity;
    stim_info.FB = FB;
    stim_info.FB_motion = FB_motion;
    stim_info.LVL1 = LVL1;
    stim_info.LVL2 = LVL2;
    stim_info.LVL3 = LVL3;
    stim_info.LVL4 = LVL4;
    
    if ~exist('select2','var')
        select2 = nan(1,length(select));
        tdays2 = nan(1,length(select));

    end
    select2(select == 0) = [];
    select2(select2 == 0) = nan;
    select(select == 0) = [];
    tdays = 1:length(select);
    
    
    smooth_select = smooth(select,3);
    smooth_MeanDT = smooth(MeanDecisionTime,3);
    
    %       if input.doMask
        tdays2 = 1:length(select);
%         notnanidx = find(~isnan(select2));
        counter = 0;
        for j = 1:length(select2)
            if  isnan(select2(j)) 
                continue
            elseif ~isnan(select2(j)) && counter <5
                counter = counter+1;             
            else
                startidx = j;
                break
            end
        end
        
    if  exist('startidx','var') %%
        select2_temp = smooth(select2(startidx:end),3);
        smooth_select2 = nan(1,length(tdays2));
        smooth_select2(startidx:end) = select2_temp; 
    else
        smooth_select2 = nan(1,length(tdays2));

    end
    
    
%     axI = subplot(spSz{:}, 8);	
%     plot(tdays, smooth_select, 'k', 'LineWidth', 1.25);
%     hold on
%     ylim([-0.1 1])
%     hline(0.9)
%     hline(0)
%     xlabel('Training Days')
%     ylabel('Selectivity')
%     title('Total Performance History (MA)')
%     
% %     if input.doMask
%         tdays2 = 1:length(select);
% %         notnanidx = find(~isnan(select2));
%         counter = 0;
%         for j = 1:length(select2)
%             if  isnan(select2(j)) 
%                 continue
%             elseif ~isnan(select2(j)) && counter <5
%                 counter = counter+1;             
%             else
%                 startidx = j;
%                 break
%             end
%         end
%         
%     if  exist('startidx','var') %%
%         select2_temp = smooth(select2(startidx:end),3);
%         smooth_select2 = nan(1,length(tdays2));
%         smooth_select2(startidx:end) = select2_temp;
%         plot(tdays2, smooth_select2, 'b', 'LineWidth', 1.25);
%     end
%    
%     yyaxis right
%     plot(tdays, smooth_MeanDT', 'LineWidth', 1.24)
%     ylabel('Decision Time')