clear all
close all

%% data concatenation
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\';
load(fullfile(data_path, 'CDMouseList.mat'))

% mouse_names_PV_contrast = [{'i419'}; 'i548'; 'i565'; 'i578'; 'i593']
% All Curve stats

quad_data = struct()



% for m =1:length(mouse_names_PV_contrast)
%     mouse = mouse_names_PV_contrast{m}
%     load(fullfile(data_path,[mouse 'goodData.mat']));
%     
%     
%    tempStruct(m).trialOutcomeCell = input.trialOutcomeCell;
%    tempStruct(m).tDecisionTimeMs = input.tDecisionTimeMs;
%    tempStruct(m).tBlock2TrialNumber = input.tBlock2TrialNumber;
%    tempStruct(m).tGratingEccentricityDeg = input.tGratingEccentricityDeg;
%    tempStruct(m).feedbackMotionSensitivity = input.feedbackMotionSensitivity;
%    tempStruct(m).tGratingContrast = input.tGratingContrast;
%    tempStruct(m).dGratingContrast  = input.dGratingContrast;
%    tempStruct(m).tLeftTrial = input.tLeftTrial;
%    tempStruct(m).tThisTrialStartTimeMs = input.tThisTrialStartTimeMs;
%    tempStruct(m).stimTimestampMs = input.stimTimestampMs;
%    tempStruct(m).tDecisionTimeMs = input.tDecisionTimeMs;
%    tempStruct(m).block2TrialLaserPowerMw = input.block2TrialLaserPowerMw;
%    tempStruct(m).quadratureTimesUs  = input.quadratureTimesUs;
%    tempStruct(m).quadratureValues  = input.quadratureValues;
%    tempStruct(m).delayTimeMs = input.delayTimeMs;
%  
    


load(fullfile(data_path,['i461' 'goodData.mat']));

% input = concatenateDataBlocks(tempStruct); % Using input from "mouse"goodData.mat - CLM

%%  
% cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
% savedate=input.saveTime; 
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell(Tix)))); %%? Why 18001 (total time window? 8000 + 10000)? Why use uint16 (changes class to unit16 - 16-bit intergers)? - CLM
% qVals_final = qVals_final(:, Tix);
% qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
% qTimes_act = qTimes_act(:, Tix);
% qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
% qTimes_thresh = qTimes_thresh(1, Tix);
SIx = double(strcmp(input.trialOutcomeCell,'success'));
FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
block2=celleqel2mat_padded(input.tBlock2TrialNumber);
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg); % 
% eccentricity = eccentricity(Tix);
feedbackMotion = input.feedbackMotionSensitivity;
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity(1); rdt=-1*ldt; % ldt/rdt = left/right decision threshold - CLM 
% tContrast=celleqel2mat_padded(input.tGratingContrast); % comment out for size discrim
% unqtargets=unique(tContrast);
% dContrast=celleqel2mat_padded(input.dGratingContrast); % comment out for size discrim
% contrastratio=tContrast./dContrast;
% contrastratio=round(contrastratio,3);
% uncontrasts=unique(contrastratio);
uncontrasts = ratios(5:8);
ratio_mat(ratio_mat == ratios(1)) = ratios(8);
ratio_mat(ratio_mat == ratios(2)) = ratios(7);
ratio_mat(ratio_mat == ratios(3)) = ratios(6);
ratio_mat(ratio_mat == ratios(4)) = ratios(5);
contrastratio = ratio_mat;
c1=uncontrasts(4); c2=uncontrasts(3); c3=uncontrasts(2); c4=uncontrasts(1);
tLeftTrial=celleqel2mat_padded(input.tLeftTrial);
delay = double(input.delayTimeMs);
% thistrialtime=celleqel2mat_padded(input.tThisTrialStartTimeMs);
% stimtime=celleqel2mat_padded(input.stimTimestampMs);
% dectime=celleqel2mat_padded(input.tDecisionTimeMs);
% ledpwr=input.block2TrialLaserPowerMw;
% quadTimes = input.quadratureTimesUs;
% quadValues = input.quadratureValues;

% currIdx = 0;

for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
%         currIdx = currIdx + 1;
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000); %  Times of each registered wheel movement starting from current and ending of next trial. Microseconds to milliseconds
        qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]); % Wheel positions from current trial to end of next. Why is this done this way?
%           qTimes = double(input.quadratureTimesUs{trN}./1000);
%           qVals = double(input.quadratureValues{trN});
        stimTime = double(input.stimTimestampMs{trN});
%         stimTime= stimtime+delay;
        qTimes_zero = qTimes-stimTime; % Wheel movement before (-) and after (+) stim appears
        qVals = qVals-qVals(1); % Start quad values from 0
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000); % Index of times where wheel movement falls in time window
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind); % Time of wheel movement relative, to stimulus appearance, within time window
            qVals_sub = qVals(time_ind); % Wheel position in time window
%             qTimes_temp = qTimes(time_ind); % Time of raw (non-aligned) wheel movement in time window
            rep_ind = find(diff(qTimes_sub)==0); % ? Returns index where the difference between adjacent elements of q_times_sub is 0. Repeated time (end and beginning of next trial?)
            qTimes_sub(rep_ind) = []; % removes rep_ind (repeated index) from vector
            qVals_sub(rep_ind) = [];
%             qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000; % Time window
%             qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)'; % Not sure why this is being done
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)'; % interpolated movement
            if input.tDecisionTimeMs{trN} < 10000
                if isnan(qVals_final(8000,trN))
                    data_idx = find(~isnan(qVals_final(:,trN)));
                    shift_idx = 8000 + length(data_idx);
                    qVals_final(8000:shift_idx-1,trN) = qVals_final(data_idx, trN); % move data to start at 8000 (stim-on)
                    qVals_final(shift_idx:end,trN) = NaN;
                end
%                 qVal_off = qVals_final(:,trN)-qVals_final(8000,trN); %align final quad values to 0, currently not used
%                 qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN); % currently not used...
            end
        else
            continue
        end
    end
end

% trN = 125
% qTimes_test = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
% qVals_test = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]);
% stimTime = double(input.stimTimestampMs{trN}) + 800; % need +time till stimulus first moves
% qTimes_zero_test = qTimes_test-stimTime;
% dectime = stimTime+double(input.tDecisionTimeMs{trN});
% idectime = find(qTimes_test < dectime, 1, 'last')
% figure;
% plot(qTimes_test, qVals_test)
% hold on
% vline(stimTime, 'k:')
% vline(dectime, 'r:')

%%
if sum(block2)==0
    nledcond=1;
else
    nledcond=2;
end
ncorrectcond=2;
nstimlocation=2;
nratios=length(uncontrasts);
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast

FIx = FIx(1:length(FIx)-1);
SIx = SIx(1:length(SIx)-1);
tLeftTrial = tLeftTrial(1:length(tLeftTrial)-1);
contrastratio = contrastratio(1:length(contrastratio)-1);

sorted_trials{nledcond, ncorrectcond, nstimlocation, nratios} = {};
for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear tempcon tempbehave temploc temprat tempIdx
               tempcon = find(block2 == icontrol-1);
               if ibehavior == 1
                   tempbehave = find(FIx == 1);
               elseif ibehavior == 2
                   tempbehave = find(SIx == 1);
               end
               temploc = find(tLeftTrial == iloc-1);
               temprat = find(contrastratio == uncontrasts(irat));
               tempIdx =intersect(intersect(intersect(tempcon, tempbehave), temploc), temprat); % indices of trials where these conditions are aligned
               sorted_trials{icontrol, ibehavior, iloc, irat} = tempIdx; % indices of trials aligned into matrix by condition
           end
       end
   end
end



% sorted_trials = nan(nledcond, ncorrectcond, nstimlocation, nratios, numtrials);
% for icontrol=1:nledcond
%    for ibehavior = 1:ncorrectcond
%        for iloc = 1:nstimlocation
%            for irat = 1:nratios
%                clear tempcon tempbehavesuccess tempbehavefailure tempbehave temploc temprat temp
%                tempcon = find(block2 == icontrol-1);
%                if ibehavior == 1
%                    tempbehave = find(FIx == 1);
%                elseif ibehavior == 2
%                    tempbehave = find(SIx == 1);
%                end
%                temploc = find(tLeftTrial == iloc-1);
%                temprat = find(contrastratio == uncontrasts(irat));
%                temp =intersect(intersect(intersect(tempcon, tempbehave), temploc), temprat); % indexes of trials where these conditions are aligned
%                sorted_trials(icontrol, ibehavior, iloc, irat, 1:length(temp)) = temp; % sorted these aligned conditons into matrix by condition
%            end; end; end; end
%%
bluecolors=flipud(brewermap(11, '*Blues'));
redcolors=flipud(brewermap(11, '*Reds'));
rxntime=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxntime_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxntime_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);

qvals_adj=qVals_final-qVals_final(8000,:);    

rc_empty_control = 0;
rc_full_control = 0;
ri_empty_control = 0;
ri_full_control = 0;
lc_empty_control = 0;
lc_full_control = 0;
li_empty_control = 0;
li_full_control = 0;

rc_empty_LED = 0;
rc_full_LED = 0;
ri_empty_LED = 0;
ri_full_LED = 0;
lc_empty_LED = 0;
lc_full_LED = 0;
li_empty_LED = 0;
li_full_LED = 0;

missing_trials = [];
control_trials = [];
LED_trials = [];

for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear temp qvals_temp_mean time_adj temptime_itrial temptime
                if icontrol==1
                       usecolor=bluecolors;
                   elseif icontrol==2
                       usecolor=redcolors; 
                end
               condIdx = sorted_trials{icontrol,ibehavior,iloc,irat};               
%                qvals_adj=qVals_final(:,condIdx)-qVals_final(8000,condIdx);               
%                qvals_temp_mean=nanmean(qvals_adj,2);  % mean qVals of trials that share same conditions
               qvals_temp_mean=nanmean(qvals_adj(:, condIdx),2); % try this
               temptime = [ ];
                   for itrial=1:length(condIdx) 
                       iqvals_trial = condIdx(itrial); %qvals adj column index (trial)
                       iqvals = qvals_adj(8000:18000,iqvals_trial); % qvals for trial
                      if iloc==1 && ibehavior==2 %right correct (neg values)
                          ifirst_qval = find(iqvals <=rdt,1,'first');
                          if isempty(ifirst_qval) % Why does this happen if this is a right correct trial...
                              if icontrol == 1
                                  rc_empty_control = rc_empty_control+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(1)
%                                   subplot(2,2,1)
%                                   hold on
%                                   title('Empty Right Correct Control Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([-500 6000])
                              else
                                  rc_empty_LED = rc_empty_LED+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(1)
%                                   subplot(2,2,2)
%                                   hold on
%                                   title('Empty Right Correct LED Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([-500 6000])
                              end
%                               
                              continue
                          else
                              if icontrol == 1
                                  rc_full_control = rc_full_control+1;
                                  control_trials = [control_trials condIdx(itrial)];
                        
                              else
                                  rc_full_LED = rc_full_LED+1;
                                  LED_trials = [LED_trials condIdx(itrial)];
                         
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end  
                      elseif iloc==1 && ibehavior==1 %right incorrect (pos values)
                          ifirst_qval = find(iqvals >=ldt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  ri_empty_control = ri_empty_control+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(1)
%                                   subplot(2,2,3)
%                                   hold on
%                                   title('Empty Right Incorrect Control Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              else
                                  ri_empty_LED = ri_empty_LED+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(1)
%                                   subplot(2,2,4)
%                                   hold on
%                                   title('Empty Right Inorrect LED Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  ri_full_control = ri_full_control+1;
                                  control_trials = [control_trials condIdx(itrial)];
                              else
                                  ri_full_LED = ri_full_LED+1;
                                  LED_trials = [LED_trials condIdx(itrial)];
                              end
                              temptime_itrial = ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end                      
                      elseif iloc==2 && ibehavior==2 %left correct (pos values)
                          ifirst_qval = find(iqvals >=ldt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  lc_empty_control = lc_empty_control+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(2)
%                                   subplot(2,2,1)
%                                   hold on
%                                   title('Empty Left Correct Control Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              else
                                  lc_empty_LED = lc_empty_LED+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(2)
%                                   subplot(2,2,2)
%                                   hold on
%                                   title('Empty Left Correct LED Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  lc_full_control = lc_full_control+1;
                                  control_trials = [control_trials condIdx(itrial)];
                              else
                                  lc_full_LED = lc_full_LED+1;
                                  LED_trials = [LED_trials condIdx(itrial)];
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end
                      elseif iloc==2 && ibehavior==1 %left incorrect (neg values)
                          ifirst_qval = find(iqvals <=rdt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  li_empty_control = li_empty_control+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(2)
%                                   subplot(2,2,3)
%                                   hold on
%                                   title('Empty Left Inorrect Control Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              else
                                  li_empty_LED = li_empty_LED+1;
                                  missing_trials = [missing_trials condIdx(itrial)];
%                                   figure(2)
%                                   subplot(2,2,4)
%                                   hold on
%                                   title('Empty Left Inorrect LED Trials')
%                                   plot(iqvals)
%                                   ylim([-100 100]); xlim([0 6000])
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  li_full_control = li_full_control+1;
                                  control_trials = [control_trials condIdx(itrial)];
                              else
                                  li_full_LED = li_full_LED+1;
                                  LED_trials = [LED_trials condIdx(itrial)];
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end   
                      end
                      
                   rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))=temptime; %times (where quad reached threshold) across trials matching conditions
                   rxntime_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(temptime); % average times (where quad reached threshold) matching conditions
                   
                   rxnslopewindowstart=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))-50;
                   rxnslopewindowend=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))+400; 
                   
                   rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime))=qvals_temp_mean(temptime+8000);    
                   rxnqval_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime)));

%                    qvals_adj_window = qvals_adj(8000:18000, :);
                   rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime))=abs((qvals_adj(rxnslopewindowend+8000)-qvals_adj(rxnslopewindowstart+8000)))/450;
                   rxnslope_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   
%                    rxnslope_succ_fail_avg(icontrol,:,iloc,irat,1)= mean(nanmean(rxnslope(icontrol,:,iloc,irat,1:length(temptime)))); % combined correct, incorrect trials
                   
                   rxntime_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxnslope_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));      
                   
                   rxntime_SEM(icontrol,ibehavior,iloc,irat) = rxntime_error(icontrol,ibehavior,iloc,irat)/sqrt((length(temptime)));
                   rxnslope_SEM(icontrol,ibehavior,iloc,irat) = rxnslope_error(icontrol,ibehavior,iloc,irat)/sqrt((length(temptime)));
                   end
                
%                if iloc==1 && ibehavior==2 && icontrol==1
%                    subplot(2,4,3)
%                    hold on
%                    title('Right Correct Control Trials')
%                    ylim([-200 100]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1); %
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==1 && ibehavior==2 && icontrol==2
%                    subplot(2,4,4)
%                    hold on
%                    title('Right Correct LED Trials')
%                    ylim([-200 100]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==2 && ibehavior==2 && icontrol==1
%                    subplot(2,4,1)
%                    hold on
%                    title('Left Correct Control Trials')
%                    ylim([-100 200]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                     line([0 0], [-500 500],'Color','k')
%                    line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==2 && ibehavior==2 && icontrol==2
%                    subplot(2,4,2)
%                    hold on
%                    title('Left Correct LED Trials')
%                    ylim([-100 300]); xlim([-500 6000]);
%                    plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==1 && ibehavior==1 && icontrol==1
%                    subplot(2,4,7)
%                    hold on
%                    title('Right Incorrect Control Trials')
%                    ylim([-100 200]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==1 && ibehavior==1 && icontrol==2
%                    subplot(2,4,8)
%                    hold on
%                    title('Right Incorrect LED Trials')
%                    ylim([-100 200]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==2 && ibehavior==1 && icontrol==1
%                    subplot(2,4,5)
%                    hold on
%                    title('Left Incorrect Control Trials')
%                    ylim([-200 100]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k')
%                    line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                elseif iloc==2 && ibehavior==1 && icontrol==2
%                    subplot(2,4,6)
%                    hold on
%                    title('Left Incorrect LED Trials')
%                    ylim([-200 100]); xlim([-500 6000]);
%                    checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
%                    if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
%                        continue
%                    else
%                        plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
%                    end
%                    line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
%                    scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
%                    end;
%                hold on
           end; end; end; end         

% percent_rc_control = rc_empty_control / (rc_empty_control + rc_full_control)
% percent_ri_control = ri_empty_control/ (ri_empty_control + ri_full_control)
% percent_lc_control = lc_empty_control/ (lc_empty_control + lc_full_control)
% percent_li_control = li_empty_control/ (li_empty_control + li_full_control)
% 
% percent_rc_LED = rc_empty_LED / (rc_empty_LED + rc_full_LED)
% percent_ri_LED = ri_empty_LED/ (ri_empty_LED + ri_full_LED)
% percent_lc_LED = lc_empty_LED/ (lc_empty_LED + lc_full_LED)
% percent_li_LED = li_empty_LED/ (li_empty_LED + li_full_LED)
% 
% percent.rc_control = percent_rc_control
% percent.ri_control = percent_ri_control
% percent.lc_control = percent_lc_control
% percent.li_control= percent_li_control
% percent.rc_LED= percent_rc_LED
% percent.ri_LED= percent_ri_LED
% percent.lc_LED= percent_lc_LED
% percent.li_LED= percent_li_LED


Names = [{'fail'}; {'pass'}; {'% fail'}]
% Rc = [rc_empty_control; rc_full_control; percent_rc_control* 100]
% Ri = [ri_empty_control; ri_full_control; percent_ri_control* 100]
% Lc = [lc_empty_control; lc_full_control; percent_lc_control* 100]
% Li = [li_empty_control; li_full_control; percent_li_control* 100]
% 
% Rc_LED = [rc_empty_LED; rc_full_LED; percent_rc_LED* 100]
% Ri_LED = [ri_empty_LED; ri_full_LED; percent_ri_LED* 100]
% Lc_LED = [lc_empty_LED; lc_full_LED; percent_lc_LED* 100]
% Li_LED = [li_empty_LED; li_full_LED; percent_li_LED* 100]
% 
% error_table_by side = table(Names, Rc, Rc_LED,  Ri, Ri_LED,  Lc, Lc_LED, Li, Li_LED)



percent_R_control = (rc_empty_control + ri_empty_control) / (rc_empty_control + rc_full_control + ri_empty_control + ri_full_control) * 100;
percent_L_control = (lc_empty_control + li_empty_control)/ (lc_empty_control + lc_full_control + li_empty_control + li_full_control) * 100;

percent_R_LED = (rc_empty_LED + ri_empty_LED) / (rc_empty_LED + rc_full_LED + ri_empty_LED + ri_full_LED) * 100;
percent_L_LED = (lc_empty_LED + li_empty_LED)/ (lc_empty_LED + lc_full_LED + li_empty_LED + li_full_LED) * 100;

R = [(rc_empty_control + ri_empty_control); (rc_full_control + ri_full_control); percent_R_control ];
R_LED = [(rc_empty_LED + ri_empty_LED); (rc_full_LED + ri_full_LED); percent_R_LED ];
L = [(lc_empty_control + li_empty_control); (lc_full_control + li_full_control); percent_L_control ];
L_LED = [(lc_empty_LED + li_empty_LED); (lc_full_LED + li_full_LED); percent_L_LED ];

error_table = table(Names, R, R_LED, L, L_LED);


%% Find %correct by condition

%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast

sub_sort = sorted_trials;

for i = 1:length(missing_trials)
    curr_trial = missing_trials(i);
    for icontrol=1:nledcond
        for ibehavior = 1:ncorrectcond
            for iloc = 1:nstimlocation
                for irat = 1:nratios
                    if find(sorted_trials{icontrol,ibehavior,iloc,irat} == curr_trial)
                      strip_idx = sorted_trials{icontrol,ibehavior,iloc,irat} == curr_trial;
                      trials = sorted_trials{icontrol,ibehavior,iloc,irat} ;
                      trials(strip_idx) = [];
                      sub_sort{icontrol,ibehavior,iloc,irat} = trials;
                    end
                end
            end
        end
    end
end


% No removal of trials that fail to reach threshold 
% for icontrol=1:nledcond
%    for iloc = 1:nstimlocation
%        for irat = 1:nratios
%            correct = length(sorted_trials{icontrol,2,iloc,irat});
%            incorrect = length(sorted_trials{icontrol,1,iloc,irat});
%            total = correct + incorrect;
%            prct_corr(icontrol, iloc, irat) = correct / total;
%        end
%    end
% end

%Removal of trials

for icontrol=1:nledcond
   for iloc = 1:nstimlocation
       for irat = 1:nratios
           correct(icontrol, iloc, irat) = length(sub_sort{icontrol,2,iloc,irat});
           incorrect(icontrol, iloc, irat) = length(sub_sort{icontrol,1,iloc,irat});
           total = correct + incorrect;
           prct_corr(icontrol, iloc, irat) = correct(icontrol, iloc, irat) / total(icontrol, iloc, irat);
       end
   end
end


%% 
% % scatterplot %correct (left/right) vs. slope (by contrast ratio) (All trials, Control vs LED
% % same plot
% 
% % scatterrplot %incorrect (left/right) vs. slope (by contrast ratio) (Incorrect trials, Control vs LED
% % same plot
% 

taskstr='Contrast';

%rxnslope_avg(icontrol,ibehavior,iloc,irat)
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast

for i = 1:4
    mscC_left(i) = rxnslope_avg(1,2,2,i);
    mscC_right(i) = rxnslope_avg(1,2,1,i);
    
    mscLED_left(i) = rxnslope_avg(2,2,2,i);
    mscLED_right(i) = rxnslope_avg(2,2,1,i);
    
    msiC_left(i) = rxnslope_avg(1,1,2,i);
    msiC_right(i) = rxnslope_avg(1,1,1,i);
    
    msiLED_left(i) = rxnslope_avg(2,1,2,i);
    msiLED_right(i) = rxnslope_avg(2,1,1,i);
    
    %error (std)
    
    mscC_left_error(i) = rxnslope_error(1,2,2,i);
    mscC_right_error(i) = rxnslope_error(1,2,1,i);
    
    mscLED_left_error(i) = rxnslope_error(2,2,2,i);
    mscLED_right_error(i) = rxnslope_error(2,2,1,i);
    
    msiC_left_error(i) = rxnslope_error(1,1,2,i);
    msiC_right_error(i) = rxnslope_error(1,1,1,i);
    
    msiLED_left_error(i) = rxnslope_error(2,1,2,i);
    msiLED_right_error(i) = rxnslope_error(2,1,1,i);
    
    %error (SEM)
    
    mscC_left_SEM(i) = rxnslope_SEM(1,2,2,i);
    mscC_right_SEM(i) = rxnslope_SEM(1,2,1,i);
    
    mscLED_left_SEM(i) = rxnslope_SEM(2,2,2,i);
    mscLED_right_SEM(i) = rxnslope_SEM(2,2,1,i);
    
    msiC_left_SEM(i) = rxnslope_SEM(1,1,2,i);
    msiC_right_SEM(i) = rxnslope_SEM(1,1,1,i);
    
    msiLED_left_SEM(i) = rxnslope_SEM(2,1,2,i);
    msiLED_right_SEM(i) = rxnslope_SEM(2,1,1,i);
    
    %time to thresh
    
    mtcC_left(i) = rxntime_avg(1,2,2,i);
    mtcC_right(i) = rxntime_avg(1,2,1,i);
    
    mtcLED_left(i) = rxntime_avg(2,2,2,i);
    mtcLED_right(i) = rxntime_avg(2,2,1,i);
    
    mtiC_left(i) = rxntime_avg(1,1,2,i);
    mtiC_right(i) = rxntime_avg(1,1,1,i);
    
    mtiLED_left(i) = rxntime_avg(2,1,2,i);
    mtiLED_right(i) = rxntime_avg(2,1,1,i);
    
    %error SEM
    
    mtcC_left_SEM(i) = rxntime_SEM(1,2,2,i);
    mtcC_right_SEM(i) = rxntime_SEM(1,2,1,i);
    
    mtcLED_left_SEM(i) = rxntime_SEM(2,2,2,i);
    mtcLED_right_SEM(i) = rxntime_SEM(2,2,1,i);
    
    mtiC_left_SEM(i) = rxntime_SEM(1,1,2,i);
    mtiC_right_SEM(i) = rxntime_SEM(1,1,1,i);
    
    mtiLED_left_SEM(i) = rxntime_SEM(2,1,2,i);
    mtiLED_right_SEM(i) = rxntime_SEM(2,1,1,i);
    
    %combined correct/incorrect slope
    
%     meanslope_all__control_left(i) = rxnslope_succ_fail_avg(1,1,2,i);
%     meanslope_all__control_right(i) = rxnslope_succ_fail_avg(1,1,1,i);
%     
%     meanslope_all__LED_left(i) = rxnslope_succ_fail_avg(2,1,2,i);
%     meanslope_all__LED_right(i) = rxnslope_succ_fail_avg(2,1,1,i);
    
    %prct correct
    
    prct_corr_left_control(i) = prct_corr(1,2,i);
    prct_corr_right_control(i) = prct_corr(1,1,i);
    prct_corr_left_LED(i) = prct_corr(2,2,i);
    prct_corr_right_LED(i) = prct_corr(2,1,i);
end

% %slope
% meanslope_correct_control=horzcat(fliplr(mscC_left),mscC_right);
% meanslope_correct_LED=horzcat(fliplr(mscLED_left),mscLED_right);
% 
% meanslope_incorrect_control=horzcat(fliplr(msiC_left),msiC_right);
% meanslope_incorrect_LED=horzcat(fliplr(msiLED_left),msiLED_right);
% 
% % meanslope_all_control = horzcat(fliplr(meanslope_all__control_left),meanslope_all__control_right);
% % meanslope_all_LED = horzcat(fliplr(meanslope_all__LED_left),meanslope_all__LED_right);
% 
% %error
% 
% meanslope_correct_control_error=horzcat(fliplr(mscC_left_error),mscC_right_error);
% meanslope_correct_LED_error=horzcat(fliplr(mscLED_left_error),mscLED_right_error);
% 
% meanslope_incorrect_control_error=horzcat(fliplr(msiC_left_error),msiC_right_error);
% meanslope_incorrect_LED_error=horzcat(fliplr(msiLED_left_error),msiLED_right_error);
% 
% 
% 
% %time
% 
% meantime_correct_control=horzcat(fliplr(mtcC_left),mtcC_right);
% meantime_correct_LED=horzcat(fliplr(mtcLED_left),mtcLED_right);
% 
% meantime_incorrect_control=horzcat(fliplr(mtiC_left),mtiC_right);
% meantime_incorrect_LED=horzcat(fliplr(mtiLED_left),mtiLED_right);

slope_control = mean([mscC_left; mscC_right]);
slope_LED = mean([mscLED_left; mscLED_right]);

slope_control_SEM = mean([mscC_left_SEM; mscC_right_SEM]);
slope_LED_SEM = mean([mscLED_left_SEM; mscLED_right_SEM]);


time_control = mean([mtcC_left; mtcC_right]);
time_LED = mean([mtcLED_left; mtcLED_right]);

time_control_SEM = mean([mtcC_left_SEM; mtcC_right_SEM]);
time_LED_SEM = mean([mtcLED_left_SEM; mtcLED_right_SEM]);

prct_corr_control = mean([prct_corr_left_control; prct_corr_right_control]);
prct_corr_LED = mean([prct_corr_left_LED; prct_corr_right_LED]);


xright=[c4 c3 c2 c1]; xleft=[1/c1 1/c2 1/c3 1/c4];
x=horzcat(xleft,xright);
xstr=sprintf('%s ratio (R/L)',taskstr);




% figure;
% subplot(1,2,1)
% hold on
errorbar(prct_corr_left_control, mscC_left, mscC_left_error, '<k')
% errorbar(prct_corr_left_LED, mscLED_left, mscLED_left_error, '<b')
% xlim([0 1])
% xlabel("% correct")
% ylabel('slope')
% title('Left')
% 
% subplot(1,2,2)
% hold on
% errorbar(prct_corr_right_control, mscC_right, mscC_right_error, '>k')
% errorbar(prct_corr_right_LED, mscLED_right, mscLED_right_error, '>b')
% xlim([0 1])
% xlabel("% correct")
% ylabel('slope')
% title('Right')



% figure;
% for i = 1:4
%     subplot(1,4,i)
%     title(num2str(uncontrasts(i)))
%     xlabel("% correct")
%     ylabel('slope')
%     hold on
%     errorbar(prct_corr_left_control(i), mscC_left(i), mscC_left_error(i), '<k')
%     errorbar(prct_corr_left_LED(i), mscLED_left(i), mscLED_left_error(i), '<b')
%     errorbar(prct_corr_right_control(i), mscC_right(i), mscC_right_error(i), '>k')
%     errorbar(prct_corr_right_LED(i), mscLED_right(i), mscLED_right_error(i), '>b')
%     xlim([0 1])
%     ylim([-0.05 0.3])
% end


figure;
subplot(2,2,1)
% for i = 1:4
%     subplot(1,4,i)
    title([mouse ' Slope Left & Right'])
    xlabel("% correct")
    ylabel('slope')
    hold on
    errorbar(prct_corr_left_control, mscC_left, mscC_left_SEM, '<k')
    errorbar(prct_corr_left_LED, mscLED_left, mscLED_left_SEM, '<b')
    errorbar(prct_corr_right_control, mscC_right, mscC_right_SEM, '>k')
    errorbar(prct_corr_right_LED, mscLED_right, mscLED_right_SEM, '>b')
%     vline(0.5)
%     xlim([0 1])
%     ylim([-0.05 0.3])
% end

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_slope_by_correct_SEM_all_ratios'], '-dpdf')
% savefig([mouse '_slope_by_correct_SEM_all_ratios'])

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_slope_by_correct_SEM'], '-dpdf')
% savefig([mouse '_slope_by_correct_SEM'])


% figure;
subplot(2,2,2)

% for i = 1:4
%     subplot(1,4,i)
    title([mouse ' Reaction Time Left & Right'])
    xlabel("% correct")
    ylabel('RT')
    hold on
    errorbar(prct_corr_left_control, mtcC_left, mtcC_left_SEM, '<k')
    errorbar(prct_corr_left_LED, mtcLED_left, mtcLED_left_SEM, '<b')
    errorbar(prct_corr_right_control, mtcC_right, mtcC_right_SEM, '>k')
    errorbar(prct_corr_right_LED, mtcLED_right, mtcLED_right_SEM, '>b')
%     vline(0.5)
%     xlim([0 1])
%     ylim([-0.05 0.3])
% end

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_RT_by_correct_SEM_all_ratios'], '-dpdf')
% savefig([mouse '_RT_by_correct_SEM_all_ratios'])


% figure;
subplot(2,2,3)

% for i = 1:4
%     subplot(1,4,i)
    title([mouse ' Slope'])
    xlabel("% correct")
    ylabel('slope')
    hold on
    errorbar(prct_corr_control, slope_control, slope_control_SEM, 'ok')
    errorbar(prct_corr_LED, slope_LED, slope_LED_SEM, 'ob')
%     vline(0.5)
%     xlim([0 1])
%     ylim([-0.05 0.3])
% end

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_slope_by_correct'], '-dpdf')
% savefig([mouse '_slope_by_correct'])

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_slope_by_correct_SEM'], '-dpdf')
% savefig([mouse '_slope_by_correct_SEM'])


% figure;
subplot(2,2,4)

% for i = 1:4
%     subplot(1,4,i)
    title([mouse ' Reaction Time'])
    xlabel("% correct")
    ylabel('RT')
    hold on
    errorbar(prct_corr_control, time_control, time_control_SEM, 'ok')
    errorbar(prct_corr_LED, time_LED, time_LED_SEM, 'ob')

%     vline(0.5)
%     xlim([0 1])
%     ylim([-0.05 0.3])
% end

cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
print([mouse '_Slope_RT_by_correct_ratio'], '-dpdf')
% savefig([mouse '_RT_by_correct'])



quad_data(m).mouse = mouse;
quad_data(m).uncontrasts = uncontrasts;
quad_data(m).rxntime = rxntime;
quad_data(m).rxntime_avg = rxntime_avg;
quad_data(m).rxnqval = rxnqval;
quad_data(m).rxnqval_avg = rxnqval_avg;
quad_data(m).rxnslope = rxnslope;
quad_data(m).rxnslope_avg = rxnslope_avg;
quad_data(m).rxntime_error = rxntime_error;
quad_data(m).rxnslope_error = rxnslope_error;
quad_data(m).error_table = error_table;
quad_data(m).sorted_trials = sorted_trials;
quad_data(m).sub_sort = sub_sort;
quad_data(m).missing_trials = missing_trials;
quad_data(m).prct_corr_left_control = prct_corr_left_control;
quad_data(m).mscC_left = mscC_left;
quad_data(m).mscC_left_error = mscC_left_error;
quad_data(m).prct_corr_left_LED = prct_corr_left_LED;
quad_data(m).mscLED_left = mscLED_left;
quad_data(m).mscLED_left_error =  mscLED_left_error;
quad_data(m).prct_corr_right_control = prct_corr_right_control;
quad_data(m).mscC_right = mscC_right;
quad_data(m).mscC_right_error = mscC_right_error;
quad_data(m).prct_corr_right_LED = prct_corr_right_LED;
quad_data(m).mscLED_right = mscLED_right;
quad_data(m).mscLED_right_error = mscLED_right_error;
quad_data(m).prct_corr = prct_corr;
quad_data(m).correct = correct;
quad_data(m).incorrect = incorrect;
quad_data(m).rxnslope_SEM = rxnslope_SEM;
quad_data(m).mscC_left_SEM = mscC_left_SEM;
quad_data(m).mscC_right_SEM = mscC_right_SEM;
quad_data(m).mscLED_left_SEM = mscLED_left_SEM;
quad_data(m).mscLED_right_SEM = mscLED_right_SEM;


% figure;
% scatter(x ,meanslope_correct_control, 'k')
% hold on
% scatter(x ,meanslope_correct_LED, 'b')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]);
% set(gca,'box','off')
% xlabel(xstr), ylabel('Slope'); 
% 

% 
% figure;
% % figure(2)
% subplot(2,2,3)
% scatter(x, meanslope_correct_control, 'k')
% % errorbar(x, meanslope_correct_control, meanslope_correct_control_error,'ok')
% hold on
% scatter(x, meanslope_correct_LED, 'b')
% % errorbar(x, meanslope_correct_LED, meanslope_correct_LED_error,'ob')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]); ylim([0 0.4]);
% plot([1,1],[0,1], 'k', 'LineStyle', '--');
% xlabel(xstr); ylabel('Quadrature trace slope'); 
% titlestr=sprintf('Mean quad slope by contrast ratio - correct');
% title(titlestr)
% % legend('Control', 'LED')
% 
% 
% % figure(3)
% subplot(2,2,4)
% scatter(x, meanslope_incorrect_control, 'k')
% hold on
% scatter(x, meanslope_incorrect_LED, 'b')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]); ylim([0 0.4]);
% plot([1,1],[0,1], 'k', 'LineStyle', '--');
% xlabel(xstr); ylabel('Quadrature trace slope'); 
% titlestr=sprintf('Mean quad slope by contrast ratio - incorrect');
% title(titlestr)
% 
% % figure(4)
% subplot(2,2,1)
% scatter(x, meantime_correct_control, 'k')
% hold on
% scatter(x, meantime_correct_LED, 'b')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]); ylim([0 1000]);
% plot([1,1],[0,10000], 'k', 'LineStyle', '--');
% xlabel(xstr); ylabel('Response Time'); 
% titlestr=sprintf('Mean response time by contrast ratio - correct');
% title(titlestr)
% % legend('Control', 'LED')
% 
% 
% % figure(5)
% subplot(2,2,2)
% scatter(x, meantime_incorrect_control, 'k')
% hold on
% scatter(x, meantime_incorrect_LED, 'b')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]); ylim([0 5000]);
% plot([1,1],[0,10000], 'k', 'LineStyle', '--');
% xlabel(xstr); ylabel('Response Time'); 
% titlestr=sprintf('Mean response time by contrast ratio - incorrect');
% title(titlestr)

% figure(6)
% scatter3(x,pct_Corr_all(1,:),meanslope_all_control, 'k')
% hold on
% scatter3(x,pct_Corr_all(2,:),meanslope_all_LED, 'b')
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlim([0 102]);
% xlabel(xstr); ylabel('% Correct'), zlabel('Slope'); 
% title('Avg quad slope by %Correct and Contrast Ratio');
% 
% save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\'...
%      mouse '_quad_info.mat'], 'quad_data');


clearvars -except allmice data_path mouse_names_PV_contrast quad_data

% end

save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\'...
     mouse '_quad_info.mat'], 'quad_data');

%%
% clearvars -except rxntime_correct_right rxntime_correct_left rxntime_incorrect_right rxntime_incorrect_left...
%     rxnslope_correct_right rxnslope_correct_left rxnslope_incorrect_right rxnslope_incorrect_left...
%     meanslope_correct_right meanslope_correct_left meanslope_incorrect_right meanslope_incorrect_left...
%     meanRT_correct_right meanRT_correct_left meanRT_incorrect_right meanRT_incorrect_left...
%     SEMslope_correct_right SEMslope_correct_left SEMslope_incorrect_right SEMslope_incorrect_left...
%     SEMRT_correct_right SEMRT_correct_left SEMRT_incorrect_right SEMRT_incorrect_left ...
%     acttime_correct_right acttime_correct_left acttime_incorrect_right...
%     acttime_incorrect_left SEMacttime_correct_right SEMacttime_correct_left SEMacttime_incorrect_right...
%     SEMacttime_incorrect_left meanacttime_correct_left meanacttime_correct_right meanacttime_incorrect_right...
%     meanacttime_incorrect_left


%% Add target contrast to analysis
% 
% mouse_names_PV_contrast = [{'i419'}; 'i548'; 'i565'; 'i578'; 'i593']
% 
% 
% for m =1:length(mouse_names_PV_contrast)
%     mouse = mouse_names_PV_contrast{m}
%     load(fullfile(data_path,[mouse 'goodData.mat']));
%     
% 
% % load(fullfile(data_path,['i419' 'goodData.mat']));
% % 
% % mouse = 'i419';
% 
% %%  
% Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
% Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
% maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
% qVals_final = nan(18001, uint16(length(input.trialOutcomeCell(Tix)))); %%? Why 18001 (total time window? 8000 + 10000)? Why use uint16 (changes class to unit16 - 16-bit intergers)? - CLM
% SIx = double(strcmp(input.trialOutcomeCell,'success'));
% FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
% block2=celleqel2mat_padded(input.tBlock2TrialNumber);
% eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg); % 
% feedbackMotion = input.feedbackMotionSensitivity;
% ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity(1); rdt=-1*ldt; % ldt/rdt = left/right decision threshold - CLM 
% 
% %use contrast_mat and contrasts
% %use ratio_mat and ratios
% unratios = ratios(5:8);
% ratio_mat(ratio_mat == ratios(1)) = ratios(8);
% ratio_mat(ratio_mat == ratios(2)) = ratios(7);
% ratio_mat(ratio_mat == ratios(3)) = ratios(6);
% ratio_mat(ratio_mat == ratios(4)) = ratios(5);
% contrastratio = ratio_mat;
% c1=unratios(4); c2=unratios(3); c3=unratios(2); c4=unratios(1);
% tLeftTrial=celleqel2mat_padded(input.tLeftTrial);
% delay = double(input.delayTimeMs);
% 
% 
% % currIdx = 0;
% 
% for trN = 1:length(input.trialOutcomeCell)-1
%     if find(Tix == trN)
%         qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000); %  Times of each registered wheel movement starting from current and ending of next trial. Microseconds to milliseconds
%         qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]); % Wheel positions from current trial to end of next. Why is this done this way?
%         stimTime = double(input.stimTimestampMs{trN});
%         qTimes_zero = qTimes-stimTime; % Wheel movement before (-) and after (+) stim appears
%         qVals = qVals-qVals(1); % Start quad values from 0
%         time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000); % Index of times where wheel movement falls in time window
%         if length(time_ind)>2
%             qTimes_sub = qTimes_zero(time_ind); % Time of wheel movement relative, to stimulus appearance, within time window
%             qVals_sub = qVals(time_ind); % Wheel position in time window
%             rep_ind = find(diff(qTimes_sub)==0); % ? Returns index where the difference between adjacent elements of q_times_sub is 0. Repeated time (end and beginning of next trial?)
%             qTimes_sub(rep_ind) = []; % removes rep_ind (repeated index) from vector
%             qVals_sub(rep_ind) = [];
%             qTimes_final = -8000:10000; % Time window
%             qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)'; % interpolated movement
%             if input.tDecisionTimeMs{trN} < 10000
%                 if isnan(qVals_final(8000,trN))
%                     data_idx = find(~isnan(qVals_final(:,trN)));
%                     shift_idx = 8000 + length(data_idx);
%                     qVals_final(8000:shift_idx-1,trN) = qVals_final(data_idx, trN); % move data to start at 8000 (stim-on)
%                     qVals_final(shift_idx:end,trN) = NaN;
%                 end
% 
%             end
%         else
%             return
%         end
%     end
% end
% 
% 
% %%
% if sum(block2)==0
%     nledcond=1;
% else
%     nledcond=2;
% end
% ncorrectcond=2;
% nstimlocation=2;
% nratios=length(unratios);
% ncontrasts = length(contrasts);
% %control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
% %irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast
% 
% FIx = FIx(1:length(FIx)-1);
% SIx = SIx(1:length(SIx)-1);
% tLeftTrial = tLeftTrial(1:length(tLeftTrial)-1);
% contrastratio = contrastratio(1:length(contrastratio)-1);
% 
% sorted_trials{nledcond, ncorrectcond, nstimlocation, nratios, ncontrasts} = {};
% for icontrol=1:nledcond
%    for ibehavior = 1:ncorrectcond
%        for iloc = 1:nstimlocation
%            for irat = 1:nratios
%                for icon = 1:ncontrasts
%                    clear tempcon tempbehave temploc temprat tempIdx
%                    tempLED = find(block2 == icontrol-1);
%                    if ibehavior == 1
%                        tempbehave = find(FIx == 1);
%                    elseif ibehavior == 2
%                        tempbehave = find(SIx == 1);
%                    end
%                    temploc = find(tLeftTrial == iloc-1);
%                    temprat = find(contrastratio == unratios(irat));
%                    tempcon = find(contrast_mat == contrasts(icon));
%                    tempIdx =intersect(intersect(intersect(intersect(tempLED, tempbehave), temploc), temprat), tempcon); % indices of trials where these conditions are aligned
%                    sorted_trials{icontrol, ibehavior, iloc, irat, icon} = tempIdx; % indices of trials aligned into matrix by condition
%                end
%            end
%        end
%    end
% end
% 
% %%
% bluecolors=flipud(brewermap(11, '*Blues'));
% redcolors=flipud(brewermap(11, '*Reds'));
% rxntime=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxntime_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxnqval=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxnqval_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxnslope=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxnslope_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxntime_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% rxnslope_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,ncontrasts,1);
% 
% qvals_adj=qVals_final-qVals_final(8000,:);    
% 
% rc_empty_control = 0;
% rc_full_control = 0;
% ri_empty_control = 0;
% ri_full_control = 0;
% lc_empty_control = 0;
% lc_full_control = 0;
% li_empty_control = 0;
% li_full_control = 0;
% 
% rc_empty_LED = 0;
% rc_full_LED = 0;
% ri_empty_LED = 0;
% ri_full_LED = 0;
% lc_empty_LED = 0;
% lc_full_LED = 0;
% li_empty_LED = 0;
% li_full_LED = 0;
% 
% missing_trials = [];
% control_trials = [];
% LED_trials = [];
% 
% for icontrol=1:nledcond
%    for ibehavior = 1:ncorrectcond
%        for iloc = 1:nstimlocation
%            for irat = 1:nratios
%                for icon = 1:ncontrasts
%                    clear temp qvals_temp_mean time_adj temptime_itrial temptime
%                     if icontrol==1
%                            usecolor=bluecolors;
%                        elseif icontrol==2
%                            usecolor=redcolors; 
%                     end
%                    condIdx = sorted_trials{icontrol,ibehavior,iloc,irat, icon};               
%                    qvals_temp_mean=nanmean(qvals_adj(:, condIdx),2); % try this
%                    temptime = [ ];
%                        for itrial=1:length(condIdx) 
%                            iqvals_trial = condIdx(itrial); %qvals adj column index (trial)
%                            iqvals = qvals_adj(8000:18000,iqvals_trial); % qvals for trial
%                           if iloc==1 && ibehavior==2 %right correct (neg values)
%                               ifirst_qval = find(iqvals <=rdt,1,'first');
%                               if isempty(ifirst_qval) % Why does this happen if this is a right correct trial...
%                                   if icontrol == 1
%                                       rc_empty_control = rc_empty_control+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   else
%                                       rc_empty_LED = rc_empty_LED+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   end
%     %                               
%                                   continue
%                               else
%                                   if icontrol == 1
%                                       rc_full_control = rc_full_control+1;
%                                       control_trials = [control_trials condIdx(itrial)];
% 
%                                   else
%                                       rc_full_LED = rc_full_LED+1;
%                                       LED_trials = [LED_trials condIdx(itrial)];
% 
%                                   end
%                                   temptime_itrial =ifirst_qval;
%                                   temptime = [temptime temptime_itrial];
%                               end  
%                           elseif iloc==1 && ibehavior==1 %right incorrect (pos values)
%                               ifirst_qval = find(iqvals >=ldt,1,'first');
%                               if isempty(ifirst_qval)
%                                   if icontrol == 1
%                                       ri_empty_control = ri_empty_control+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   else
%                                       ri_empty_LED = ri_empty_LED+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   end
% 
%                                   continue
%                               else
%                                   if icontrol == 1
%                                       ri_full_control = ri_full_control+1;
%                                       control_trials = [control_trials condIdx(itrial)];
%                                   else
%                                       ri_full_LED = ri_full_LED+1;
%                                       LED_trials = [LED_trials condIdx(itrial)];
%                                   end
%                                   temptime_itrial = ifirst_qval;
%                                   temptime = [temptime temptime_itrial];
%                               end                      
%                           elseif iloc==2 && ibehavior==2 %left correct (pos values)
%                               ifirst_qval = find(iqvals >=ldt,1,'first');
%                               if isempty(ifirst_qval)
%                                   if icontrol == 1
%                                       lc_empty_control = lc_empty_control+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   else
%                                       lc_empty_LED = lc_empty_LED+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   end
%     %         
%                                   continue
%                               else
%                                   if icontrol == 1
%                                       lc_full_control = lc_full_control+1;
%                                       control_trials = [control_trials condIdx(itrial)];
%                                   else
%                                       lc_full_LED = lc_full_LED+1;
%                                       LED_trials = [LED_trials condIdx(itrial)];
%                                   end
%                                   temptime_itrial =ifirst_qval;
%                                   temptime = [temptime temptime_itrial];
%                               end
%                           elseif iloc==2 && ibehavior==1 %left incorrect (neg values)
%                               ifirst_qval = find(iqvals <=rdt,1,'first');
%                               if isempty(ifirst_qval)
%                                   if icontrol == 1
%                                       li_empty_control = li_empty_control+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
%     % 
%                                   else
%                                       li_empty_LED = li_empty_LED+1;
%                                       missing_trials = [missing_trials condIdx(itrial)];
% 
%                                   end
%     %                    
%                                   continue
%                               else
%                                   if icontrol == 1
%                                       li_full_control = li_full_control+1;
%                                       control_trials = [control_trials condIdx(itrial)];
%                                   else
%                                       li_full_LED = li_full_LED+1;
%                                       LED_trials = [LED_trials condIdx(itrial)];
%                                   end
%                                   temptime_itrial =ifirst_qval;
%                                   temptime = [temptime temptime_itrial];
%                               end   
%                           end
% 
%                        rxntime(icontrol,ibehavior,iloc,irat,icon,1:length(temptime))=temptime; %times (where quad reached threshold) across trials matching conditions
%                        rxntime_avg(icontrol,ibehavior,iloc,irat,icon,1)=nanmean(temptime); % average times (where quad reached threshold) matching conditions
% 
%                        rxnslopewindowstart=rxntime(icontrol,ibehavior,iloc,irat,icon,1:length(temptime))-50;
%                        rxnslopewindowend=rxntime(icontrol,ibehavior,iloc,irat,icon,1:length(temptime))+400; 
% 
%                        rxnqval(icontrol,ibehavior,iloc,irat,icon,1:length(temptime))=qvals_temp_mean(temptime+8000);    
%                        rxnqval_avg(icontrol,ibehavior,iloc,irat,icon,1)=nanmean(rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime)));
% 
%                        qvals_adj_window = qvals_adj(8000:18000, :);
%                        rxnslope(icontrol,ibehavior,iloc,irat,icon,1:length(temptime))=abs((qvals_adj(rxnslopewindowend+8000)-qvals_adj(rxnslopewindowstart+8000)))/450;
%                        rxnslope_avg(icontrol,ibehavior,iloc,irat,icon,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,icon,1:length(temptime)));
% 
%     %                    rxnslope_succ_fail_avg(icontrol,:,iloc,irat,1)= mean(nanmean(rxnslope(icontrol,:,iloc,irat,1:length(temptime)))); % combined correct, incorrect trials
% 
%                        rxntime_error(icontrol,ibehavior,iloc,irat,icon,1)=nanstd(rxntime(icontrol,ibehavior,iloc,irat,icon,1:length(temptime)));
%                        rxnslope_error(icontrol,ibehavior,iloc,irat,icon,1)=nanstd(rxnslope(icontrol,ibehavior,iloc,irat,icon,1:length(temptime)));      
% 
%                        rxntime_SEM(icontrol,ibehavior,iloc,irat,icon) = rxntime_error(icontrol,ibehavior,iloc,irat,icon)/sqrt((length(temptime)));
%                        rxnslope_SEM(icontrol,ibehavior,iloc,irat,icon) = rxnslope_error(icontrol,ibehavior,iloc,irat,icon)/sqrt((length(temptime)));
%                        end
% 
%                end
%            end
%        end
%    end
% end         
% 
% 
% Names = [{'fail'}; {'pass'}; {'% fail'}]
% 
% 
% 
% percent_R_control = (rc_empty_control + ri_empty_control) / (rc_empty_control + rc_full_control + ri_empty_control + ri_full_control) * 100;
% percent_L_control = (lc_empty_control + li_empty_control)/ (lc_empty_control + lc_full_control + li_empty_control + li_full_control) * 100;
% 
% percent_R_LED = (rc_empty_LED + ri_empty_LED) / (rc_empty_LED + rc_full_LED + ri_empty_LED + ri_full_LED) * 100;
% percent_L_LED = (lc_empty_LED + li_empty_LED)/ (lc_empty_LED + lc_full_LED + li_empty_LED + li_full_LED) * 100;
% 
% R = [(rc_empty_control + ri_empty_control); (rc_full_control + ri_full_control); percent_R_control ];
% R_LED = [(rc_empty_LED + ri_empty_LED); (rc_full_LED + ri_full_LED); percent_R_LED ];
% L = [(lc_empty_control + li_empty_control); (lc_full_control + li_full_control); percent_L_control ];
% L_LED = [(lc_empty_LED + li_empty_LED); (lc_full_LED + li_full_LED); percent_L_LED ];
% 
% error_table = table(Names, R, R_LED, L, L_LED);
% 
% 
% %% Find %correct by condition
% 
% %control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
% %irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast
% 
% sub_sort = sorted_trials;
% 
% for i = 1:length(missing_trials)
%     curr_trial = missing_trials(i);
%     for icontrol=1:nledcond
%         for ibehavior = 1:ncorrectcond
%             for iloc = 1:nstimlocation
%                 for irat = 1:nratios
%                     for icon = 1:ncontrasts
%                         if find(sorted_trials{icontrol,ibehavior,iloc,irat,icon} == curr_trial)
%                           strip_idx = sorted_trials{icontrol,ibehavior,iloc,irat,icon} == curr_trial;
%                           trials = sorted_trials{icontrol,ibehavior,iloc,irat,icon} ;
%                           trials(strip_idx) = [];
%                           sub_sort{icontrol,ibehavior,iloc,irat,icon} = trials;
%                         end
%                     end                         
%                 end
%             end
%         end
%     end
% end
% 
% %Removal of trials
% 
% for icontrol=1:nledcond
%    for iloc = 1:nstimlocation
%        for irat = 1:nratios
%            for icon = 1:ncontrasts
%                correct(icontrol, iloc, irat,icon) = length(sub_sort{icontrol,2,iloc,irat,icon});
%                incorrect(icontrol, iloc, irat,icon) = length(sub_sort{icontrol,1,iloc,irat,icon});
%                total = correct + incorrect;
%                prct_corr(icontrol, iloc, irat,icon) = correct(icontrol, iloc, irat,icon) / total(icontrol, iloc, irat,icon);
%            end
%        end
%    end
% end
% 
% 
% %% 
% taskstr='Contrast';
% 
% %rxnslope_avg(icontrol,ibehavior,iloc,irat, icontrast)
% %control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
% %irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest ratio and contrast
% 
% for i = 1:4
%     for j = 1:4
%     mscC_left(i,j) = rxnslope_avg(1,2,2,i,j);
%     mscC_right(i,j) = rxnslope_avg(1,2,1,i,j);
%     
%     mscLED_left(i,j) = rxnslope_avg(2,2,2,i,j);
%     mscLED_right(i,j) = rxnslope_avg(2,2,1,i,j);
%     
%     msiC_left(i,j) = rxnslope_avg(1,1,2,i,j);
%     msiC_right(i,j) = rxnslope_avg(1,1,1,i,j);
%     
%     msiLED_left(i,j) = rxnslope_avg(2,1,2,i,j);
%     msiLED_right(i,j) = rxnslope_avg(2,1,1,i,j);
%     
%     %error (std)
%     
%     mscC_left_error(i,j) = rxnslope_error(1,2,2,i,j);
%     mscC_right_error(i,j) = rxnslope_error(1,2,1,i,j);
%     
%     mscLED_left_error(i,j) = rxnslope_error(2,2,2,i,j);
%     mscLED_right_error(i,j) = rxnslope_error(2,2,1,i,j);
%     
%     msiC_left_error(i,j) = rxnslope_error(1,1,2,i,j);
%     msiC_right_error(i,j) = rxnslope_error(1,1,1,i,j);
%     
%     msiLED_left_error(i,j) = rxnslope_error(2,1,2,i,j);
%     msiLED_right_error(i,j) = rxnslope_error(2,1,1,i,j);
%     
%     %error (SEM)
%     
%     mscC_left_SEM(i,j) = rxnslope_SEM(1,2,2,i,j);
%     mscC_right_SEM(i,j) = rxnslope_SEM(1,2,1,i,j);
%     
%     mscLED_left_SEM(i,j) = rxnslope_SEM(2,2,2,i,j);
%     mscLED_right_SEM(i,j) = rxnslope_SEM(2,2,1,i,j);
%     
%     msiC_left_SEM(i,j) = rxnslope_SEM(1,1,2,i,j);
%     msiC_right_SEM(i,j) = rxnslope_SEM(1,1,1,i,j);
%     
%     msiLED_left_SEM(i,j) = rxnslope_SEM(2,1,2,i,j);
%     msiLED_right_SEM(i,j) = rxnslope_SEM(2,1,1,i,j);
%     
%     %time to thresh
%     
%     mtcC_left(i,j) = rxntime_avg(1,2,2,i,j);
%     mtcC_right(i,j) = rxntime_avg(1,2,1,i,j);
%     
%     mtcLED_left(i,j) = rxntime_avg(2,2,2,i,j);
%     mtcLED_right(i,j) = rxntime_avg(2,2,1,i,j);
%     
%     mtiC_left(i,j) = rxntime_avg(1,1,2,i,j);
%     mtiC_right(i,j) = rxntime_avg(1,1,1,i,j);
%     
%     mtiLED_left(i,j) = rxntime_avg(2,1,2,i,j);
%     mtiLED_right(i,j) = rxntime_avg(2,1,1,i,j);
%     
%     %error SEM
%     
%     mtcC_left_SEM(i,j) = rxntime_SEM(1,2,2,i,j);
%     mtcC_right_SEM(i,j) = rxntime_SEM(1,2,1,i,j);
%     
%     mtcLED_left_SEM(i,j) = rxntime_SEM(2,2,2,i,j);
%     mtcLED_right_SEM(i,j) = rxntime_SEM(2,2,1,i,j);
%     
%     mtiC_left_SEM(i,j) = rxntime_SEM(1,1,2,i,j);
%     mtiC_right_SEM(i,j) = rxntime_SEM(1,1,1,i,j);
%     
%     mtiLED_left_SEM(i,j) = rxntime_SEM(2,1,2,i,j);
%     mtiLED_right_SEM(i,j) = rxntime_SEM(2,1,1,i,j);
%     
%     %combined correct/incorrect slope
%     
% %     meanslope_all__control_left(i) = rxnslope_succ_fail_avg(1,1,2,i);
% %     meanslope_all__control_right(i) = rxnslope_succ_fail_avg(1,1,1,i);
% %     
% %     meanslope_all__LED_left(i) = rxnslope_succ_fail_avg(2,1,2,i);
% %     meanslope_all__LED_right(i) = rxnslope_succ_fail_avg(2,1,1,i);
%     
%     %prct correct
%     
%     prct_corr_left_control(i,j) = prct_corr(1,2,i,j);
%     prct_corr_right_control(i,j) = prct_corr(1,1,i,j);
%     prct_corr_left_LED(i,j) = prct_corr(2,2,i,j);
%     prct_corr_right_LED(i,j) = prct_corr(2,1,i,j);
%     end
% end

% slope_control = mean([mscC_left; mscC_right]);
% slope_LED = mean([mscLED_left; mscLED_right]);
% 
% slope_control_SEM = mean([mscC_left_SEM; mscC_right_SEM]);
% slope_LED_SEM = mean([mscLED_left_SEM; mscLED_right_SEM]);
% 
% 
% time_control = mean([mtcC_left; mtcC_right]);
% time_LED = mean([mtcLED_left; mtcLED_right]);
% 
% time_control_SEM = mean([mtcC_left_SEM; mtcC_right_SEM]);
% time_LED_SEM = mean([mtcLED_left_SEM; mtcLED_right_SEM]);

% mdl_slope_left = fitlm(prct_corr_left_control(:), mscC_left(:));
% mdl_slope_left_LED = fitlm(prct_corr_left_LED(:), mscLED_left(:));
% mdl_slope_right = fitlm(prct_corr_right_control(:), mscC_right(:));
% mdl_slope_right_LED = fitlm(prct_corr_right_LED(:), mscLED_right(:));
% 
% mdl_RT_left = fitlm(prct_corr_left_control(:), mtcC_left(:));
% mdl_RT_left_LED = fitlm(prct_corr_left_LED(:), mtcLED_left(:));
% mdl_RT_right = fitlm(prct_corr_right_control(:), mtcC_right(:));
% mdl_RT_right_LED = fitlm(prct_corr_right_LED(:), mtcLED_right(:));

% xL = prct_corr_left_control(:);
% xL_LED = prct_corr_left_LED(:);
% xR = prct_corr_right_control(:);
% xR_LED = prct_corr_right_LED(:);
% 
% S_L = mscC_left(:);
% S_L_LED = mscLED_left(:);
% S_R = mscC_right(:);
% S_R_LED = mscLED_right(:);
% RT_L = mtcC_left(:);
% RT_L_LED = mtcLED_left(:);
% RT_R = mtcC_right(:);
% RT_R_LED = mtcLED_right(:);
% 
% idx1 = isfinite(xL) & isfinite(S_L);
% idx2 = isfinite(xL_LED) & isfinite(S_L_LED);
% idx3 = isfinite(xR) & isfinite(S_R);
% idx4 = isfinite(xR_LED) & isfinite(S_R_LED);
% idx5 = isfinite(xL) & isfinite(RT_L);
% idx6 = isfinite(xL_LED) & isfinite(RT_L_LED);
% idx7 = isfinite(xR) & isfinite(RT_R);
% idx8 = isfinite(xR_LED) & isfinite(RT_R_LED);
% 
% 
% % Quad slope 
% mdl_slope_left = fitlm(xL(idx1), S_L(idx1), 'linear');
% mdl_slope_left_LED = fitlm(xL_LED(idx2), S_L_LED(idx2), 'linear');
% mdl_slope_right = fitlm(xR(idx3), S_R(idx3), 'linear');
% mdl_slope_right_LED = fitlm(xR_LED(idx4), S_R_LED(idx4), 'linear');
% 
% % RT
% 
% mdl_RT_left = fitlm(xL(idx5), RT_L(idx5), 'linear');
% outlier = isoutlier(mdl_RT_left.Residuals.Raw);
% mdl_RT_left2 = fitlm(xL(idx5), RT_L(idx5), 'linear', 'Exclude',(find(outlier)));
% 
% 
% mdl_RT_left_LED = fitlm(xL_LED(idx6), RT_L_LED(idx6), 'linear');
% outlier = isoutlier(mdl_RT_left_LED.Residuals.Raw);
% mdl_RT_left2_LED = fitlm(xL_LED(idx5), RT_L_LED(idx5), 'linear', 'Exclude',(find(outlier)));
% 
% mdl_RT_right = fitlm(xR(idx7), RT_R(idx7), 'linear');
% outlier = isoutlier(mdl_RT_right.Residuals.Raw);
% mdl_RT_right2 = fitlm(xR(idx5), RT_R(idx5), 'linear', 'Exclude',(find(outlier)));
% 
% mdl_RT_right_LED = fitlm(xR_LED(idx8), RT_R_LED(idx8), 'linear');
% outlier = isoutlier(mdl_RT_right_LED.Residuals.Raw);
% mdl_RT_right2_LED = fitlm(xR_LED(idx5), RT_R_LED(idx5), 'linear', 'Exclude',(find(outlier)));





%% 16 conditions target contrast by contrast ratio

% figure;
% % subplot(2,2,1)
% hold on
% scatter(prct_corr_left_control(:), mscC_left(:), '<k', 'filled')
% scatter(prct_corr_left_LED(:), mscLED_left(:), '<b', 'filled')
% plot(mdl_slope_left, 'k')
% plot(mdl_slope_left_LED, 'b')
% title([mouse ' Left Slope'])
% xlabel("% correct")
% ylabel('slope')
% legend('off')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_LSlope_16_points'], '-dpdf')   
%     
% 
% figure;
% % subplot(2,2,2)
% hold on
% scatter(prct_corr_right_control(:), mscC_right(:), '>k', 'filled')
% scatter(prct_corr_right_LED(:), mscLED_right(:), '>b', 'filled')
% plot(mdl_slope_right, 'k')
% plot(mdl_slope_right_LED, 'b')
% title([mouse ' Right Slope'])
% xlabel("% correct")
% ylabel('slope')
% legend('off')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_RSlope_16_points'], '-dpdf')   
% 
% 
% 
%     % time 
% 
% figure;
% % subplot(2,2,3)
% hold on
% scatter(prct_corr_left_control(:), mtcC_left(:), '<k', 'filled')
% scatter(prct_corr_left_LED(:), mtcLED_left(:), '<b', 'filled')
% plot(mdl_RT_left, 'k')
% plot(mdl_RT_left_LED, 'b')
% title([mouse ' Left RT'])
% xlabel("% correct")
% ylabel('RT')
% legend('off')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_LRT_16_points'], '-dpdf')   
% 
%     
% 
% figure;
% % subplot(2,2,4)
% hold on
% scatter(prct_corr_right_control(:), mtcC_right(:), '>k', 'filled')
% scatter(prct_corr_right_LED(:), mtcLED_right(:), '>b', 'filled')
% plot(mdl_RT_right, 'k')
% plot(mdl_RT_right_LED, 'b')
% title([mouse ' Right RT'])
% xlabel("% correct")
% ylabel('RT')
% legend('off')
% 
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_RRT_16_points'], '-dpdf')              







%%
%                
% figure;
% h1 = histogram(mscC_left, 'Normalization', 'pdf', 'FaceColor', 'k')
% hold on
% h2 = histogram(mscLED_left, 'Normalization', 'pdf', 'FaceColor', 'b')
% h1.BinWidth = 0.01;
% h2.BinWidth = 0.01;
% xlabel('slope')
% ylabel('nConditions')
% title([mouse ' Left Slope'])
% 
% y = min([h1.BinLimits h2.BinLimits]):0.001:max([h1.BinLimits h2.BinLimits]);
% mu = nanmean(mscC_left(:));
% sigma = std(mscC_left(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f, 'k', 'LineWidth',1.5)
% 
% mu = nanmean(mscLED_left(:));
% sigma = std(mscLED_left(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'b','LineWidth',1.5)
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_LSlope_hist'], '-dpdf') 
% 
% figure;
% h3= histogram(mscC_right, 'Normalization', 'pdf', 'FaceColor', 'k')
% hold on
% h4= histogram(mscLED_right, 'Normalization', 'pdf', 'FaceColor', 'b')
% h3.BinWidth = 0.01;
% h4.BinWidth = 0.01;
% xlabel('slope')
% ylabel('nConditions')
% title([mouse ' Right Slope'])
% 
% y = min([h3.BinLimits h4.BinLimits]):0.001:max([h3.BinLimits h4.BinLimits]);
% mu = nanmean(mscC_right(:));
% sigma = std(mscC_right(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f, 'k', 'LineWidth',1.5)
% 
% mu = nanmean(mscLED_right(:));
% sigma = std(mscLED_right(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'b','LineWidth',1.5)
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_RSlope_hist'], '-dpdf') 
% 
% figure;
% h5= histogram(mtcC_left, 'Normalization', 'pdf','FaceColor', 'k')
% hold on
% h6 = histogram(mtcLED_left, 'Normalization', 'pdf', 'FaceColor', 'b')
% h5.BinWidth = 100;
% h6.BinWidth = 100;
% xlabel('RT')
% ylabel('nConditions')
% title([mouse ' Left RT'])
% 
% y = min([h5.BinLimits h6.BinLimits]):0.001:max([h5.BinLimits h6.BinLimits]);
% mu = nanmean(mtcC_left(:));
% sigma = std(mtcC_left(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f, 'k', 'LineWidth',1.5)
% 
% mu = nanmean(mtcLED_left(:));
% sigma = std(mtcLED_left(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'b','LineWidth',1.5)
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_LRT_hist'], '-dpdf') 
% 
% figure;
% h7= histogram(mtcC_right, 'Normalization', 'pdf', 'FaceColor', 'k')
% hold on
% h8 = histogram(mtcLED_right, 'Normalization', 'pdf', 'FaceColor', 'b')
% h7.BinWidth = 100;
% h8.BinWidth = 100;
% xlabel('RT')
% ylabel('nConditions')
% title([mouse ' Right RT'])
% 
% y = min([h7.BinLimits h8.BinLimits]):0.001:max([h7.BinLimits h8.BinLimits]);
% mu = nanmean(mtcC_right(:));
% sigma = std(mtcC_right(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f, 'k', 'LineWidth',1.5)
% 
% mu = nanmean(mtcLED_right(:));
% sigma = std(mtcLED_right(:));
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'b','LineWidth',1.5)
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print([mouse '_RRT_hist'], '-dpdf') 

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\test\')
% save([mouse '_quadVars'], 'mdl_RT_left', 'mdl_RT_left_LED', 'mdl_RT_right', 'mdl_RT_right_LED', 'mdl_slope_left',...
% 'mdl_slope_left_LED', 'mdl_slope_right', 'mdl_slope_right_LED', 'prct_corr_left_control', 'prct_corr_left_LED', 'prct_corr_right_control',...
% 'prct_corr_right_LED', 'mscC_left', 'mscLED_left',...
% 'mscC_right','mscLED_right', 'mtcC_left', 'mtcLED_left', 'mtcC_right', 'mtcLED_right', ...
% 'mdl_RT_left2', 'mdl_RT_left2_LED', 'mdl_RT_right2', 'mdl_RT_right2_LED')
% clearvars -except allmice data_path mouse_names_PV_contrast

               
% end

%%
% data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\';
% load(fullfile(data_path, 'CDMouseList.mat'))
% mouse_names_PV_contrast = [{'i419'}; 'i548'; 'i565'; 'i578'; 'i593']
% data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\test\';
% 
% 
% 
% 
% 
% for m =1:length(mouse_names_PV_contrast)
%     mouse = mouse_names_PV_contrast{m}
%     load(fullfile(data_path,[mouse '_quadVars.mat']));
%     
% 
% p_intercept_L(m) = mdl_RT_left.Coefficients{1,4};
% p_slope_L(m) = mdl_RT_left.Coefficients{2,4};
% Rsquared_L(m) = mdl_RT_left.Rsquared.Ordinary;
% p_intercept_L2(m) = mdl_RT_left2.Coefficients{1,4};
% p_slope_L2(m) = mdl_RT_left2.Coefficients{2,4};
% Rsquared_L2(m) = mdl_RT_left2.Rsquared.Ordinary;
% 
% p_intercept_L_LED(m) = mdl_RT_left_LED.Coefficients{1,4};
% p_slope_L_LED(m) = mdl_RT_left_LED.Coefficients{2,4};
% Rsquared_L_LED(m) = mdl_RT_left_LED.Rsquared.Ordinary;
% p_intercept_L2_LED(m) = mdl_RT_left2_LED.Coefficients{1,4};
% p_slope_L2_LED(m) = mdl_RT_left2_LED.Coefficients{2,4};
% Rsquared_L2_LED(m) = mdl_RT_left2_LED.Rsquared.Ordinary;
% 
% 
% p_intercept_R(m) = mdl_RT_right.Coefficients{1,4};
% p_slope_R(m) = mdl_RT_right.Coefficients{2,4};
% Rsquared_R(m) = mdl_RT_right.Rsquared.Ordinary;
% p_intercept_R2(m) = mdl_RT_right2.Coefficients{1,4};
% p_slope_R2(m) = mdl_RT_right2.Coefficients{2,4};
% Rsquared_R2(m) = mdl_RT_right2.Rsquared.Ordinary;
% 
% p_intercept_R_LED(m) = mdl_RT_right_LED.Coefficients{1,4};
% p_slope_R_LED(m) = mdl_RT_right_LED.Coefficients{2,4};
% Rsquared_R_LED(m) = mdl_RT_right_LED.Rsquared.Ordinary;
% p_intercept_R2_LED(m) = mdl_RT_right2_LED.Coefficients{1,4};
% p_slope_R2_LED(m) = mdl_RT_right2_LED.Coefficients{2,4};
% Rsquared_R2_LED(m) = mdl_RT_right2_LED.Rsquared.Ordinary;
% 
% 
% 
% 
% slope_fitslopes_left(m) = mdl_slope_left.Coefficients{2,1};
% slope_fitslopes_leftLED(m) = mdl_slope_left_LED.Coefficients{2,1};
% slope_fitslopes_right(m) = mdl_slope_right.Coefficients{2,1};
% slope_fitslopes_rightLED(m) = mdl_slope_right_LED.Coefficients{2,1};
% 
% slope_intercepts_left(m) = mdl_slope_left.Coefficients{1,1};
% slope_intercepts_leftLED(m) = mdl_slope_left_LED.Coefficients{1,1};
% slope_intercepts_right(m) = mdl_slope_right.Coefficients{1,1};
% slope_intercepts_rightLED(m) = mdl_slope_right_LED.Coefficients{1,1};
% 
% RT_fitslopes_left(m) = mdl_RT_left.Coefficients{2,1};
% RT_fitslopes_leftLED(m) = mdl_RT_left_LED.Coefficients{2,1};
% RT_fitslopes_right(m) = mdl_RT_right.Coefficients{2,1};
% RT_fitslopes_rightLED(m) = mdl_RT_right_LED.Coefficients{2,1};
% 
% RT_intercepts_left(m) = mdl_RT_left.Coefficients{1,1};
% RT_intercepts_leftLED(m) = mdl_RT_left_LED.Coefficients{1,1};
% RT_intercepts_right(m) = mdl_RT_right.Coefficients{1,1};
% RT_intercepts_rightLED(m) = mdl_RT_right_LED.Coefficients{1,1};


% slopes_left(m) = nanmean(mscC_left(:));
% slopes_leftLED(m) = nanmean(mscLED_left(:));
% slopes_right(m) = nanmean(mscC_right(:));
% slopes_rightLED(m) = nanmean(mscLED_right(:));
% 
% RTs_left(m) = nanmean(mtcC_left(:));
% RTs_leftLED(m) = nanmean(mtcLED_left(:));
% RTs_right(m) = nanmean(mtcC_right(:));
% RTs_rightLED(m) = nanmean(mtcLED_right(:));


% end

% num_p_intercept_L = length(find(p_intercept_L <=0.05));
% num_p_slope_L = length(find(p_slope_L <=0.05));
% num_p_intercept_L2 = length(find(p_intercept_L2 <=0.05));
% num_p_slope_L2 = length(find(p_slope_L2 <=0.05));
% 
% num_p_intercept_L_LED = length(find(p_intercept_L_LED <=0.05));
% num_p_slope_L_LED = length(find(p_slope_L_LED <=0.05));
% num_p_intercept_L2_LED = length(find(p_intercept_L2_LED <=0.05));
% num_p_slope_L2_LED = length(find(p_slope_L2_LED <=0.05));
% 
% num_p_intercept_R = length(find(p_intercept_R <=0.05));
% num_p_slope_R = length(find(p_slope_R <=0.05));
% num_p_intercept_R2 = length(find(p_intercept_R2 <=0.05));
% num_p_slope_R2 = length(find(p_slope_R2 <=0.05));
% 
% num_p_intercept_R_LED = length(find(p_intercept_R_LED <=0.05));
% num_p_slope_R_LED = length(find(p_slope_R_LED <=0.05));
% num_p_intercept_R2_LED = length(find(p_intercept_R2_LED <=0.05));
% num_p_slope_R2_LED = length(find(p_slope_R2_LED <=0.05));


% figure;
% hold on
% x = categorical({'Control', 'LED'});
% x = reordercats(x,{'Control', 'LED'});
% subplot(2,2,1)
% bar(x,[[num_p_intercept_L num_p_intercept_L2]; [num_p_intercept_L_LED num_p_intercept_L2_LED]])
% ylabel('n p<=0.05')
% ylim([0 5])
% title('left intercept')
% subplot(2,2,2)
% bar(x,[[num_p_intercept_R num_p_intercept_R2]; [num_p_intercept_R_LED num_p_intercept_R2_LED]])
% ylim([0 5])
% title('right intercept')
% subplot(2,2,3)
% bar(x,[[num_p_slope_L num_p_slope_L2]; [num_p_slope_L_LED num_p_slope_L2_LED]])
% ylim([0 5])
% ylabel('n p<=0.05')
% title('left slope')
% subplot(2,2,4)
% bar(x,[[num_p_slope_R num_p_slope_R2]; [num_p_slope_R_LED num_p_slope_R2_LED]])
% ylim([0 5])
% title('right slope')
% 
% 
% 
% mean_slope_fitslopes_left = nanmean(slope_fitslopes_left);
% mean_slope_fitslopes_leftLED = nanmean(slope_fitslopes_leftLED);
% mean_slope_fitslopes_right = nanmean(slope_fitslopes_right);
% mean_slope_fitslopes_rightLED = nanmean(slope_fitslopes_rightLED);
% 
% mean_slope_intercepts_left = nanmean(slope_intercepts_left);
% mean_slope_intercepts_leftLED = nanmean(slope_intercepts_leftLED);
% mean_slope_intercepts_right = nanmean(slope_intercepts_right);
% mean_slope_intercepts_rightLED = nanmean(slope_intercepts_rightLED);
% 
% mean_RT_fitslopes_left = nanmean(RT_fitslopes_left);
% mean_RT_fitslopes_leftLED = nanmean(RT_fitslopes_leftLED);
% mean_RT_fitslopes_right = nanmean(RT_fitslopes_right);
% mean_RT_fitslopes_rightLED = nanmean(RT_fitslopes_rightLED);
% 
% mean_RT_intercepts_left = nanmean(RT_intercepts_left);
% mean_RT_intercepts_leftLED = nanmean(RT_intercepts_leftLED);
% mean_RT_intercepts_right = nanmean(RT_intercepts_right);
% mean_RT_intercepts_rightLED = nanmean(RT_intercepts_rightLED);

%%
% figure;
% hold on 
% scatter(RT_fitslopes_left, RT_fitslopes_leftLED, '<')
% scatter(RT_fitslopes_right, RT_fitslopes_rightLED, '>')
% scatter(mean_RT_fitslopes_left, mean_RT_fitslopes_leftLED, '<', 'filled')
% scatter(mean_RT_fitslopes_right, mean_RT_fitslopes_rightLED, '>', 'filled')
% xlim([-1500 400])
% ylim([-1500 400])
% refline(1, 0)
% xlabel('Control slope')
% ylabel('LED slope')
% title('RT')
% 
% figure;
% hold on
% scatter(RT_intercepts_left, RT_intercepts_leftLED, '<')
% scatter(mean_RT_intercepts_left, mean_RT_intercepts_leftLED, '<', 'filled')
% scatter(RT_intercepts_right, RT_intercepts_rightLED, '>')
% scatter(mean_RT_intercepts_right, mean_RT_intercepts_rightLED, '>', 'filled')
% xlim([200 2200])
% ylim([200 2200])
% refline(1, 0)
% xlabel('Control intercept')
% ylabel('LED intercept')
% title('RT')




% 
%%
% figure;
% hold on
% x = ones(1, length(slope_fitslopes_left));
% scatter(x*0.2, slope_fitslopes_left, '<k')
% scatter(0.2, mean(slope_fitslopes_left), '<k', 'filled')
% scatter(x*0.3, slope_fitslopes_leftLED, '<b')
% scatter(0.3, mean(slope_fitslopes_leftLED), '<b', 'filled')
% scatter(x*0.7, slope_fitslopes_right, '>k')
% scatter(0.7, mean(slope_fitslopes_right), '>k', 'filled')
% scatter(x*0.8, slope_fitslopes_rightLED, '>b')
% scatter(0.8, mean(slope_fitslopes_rightLED), '>b', 'filled')
% xlim([0 1])
% title('Slope - fit slopes')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('slope')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('slope_fitslopes', '-dpdf') 
% 
% figure;
% hold on
% x = ones(1, length(slope_intercepts_left));
% scatter(x*0.2, slope_intercepts_left, '<k')
% scatter(0.2, mean(slope_intercepts_left), '<k', 'filled')
% scatter(x*0.3, slope_intercepts_leftLED, '<b')
% scatter(0.3, mean(slope_intercepts_leftLED), '<b', 'filled')
% scatter(x*0.7, slope_intercepts_right, '>k')
% scatter(0.7, mean(slope_intercepts_right), '>k', 'filled')
% scatter(x*0.8, slope_intercepts_rightLED, '>b')
% scatter(0.8, mean(slope_intercepts_rightLED), '>b', 'filled')
% xlim([0 1])
% title('Slope - fit intercepts')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('slope')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('slope_intercepts', '-dpdf') 
%  
% figure;
% hold on
% x = ones(1, length(RT_fitslopes_left));
% scatter(x*0.2, RT_fitslopes_left, '<k')
% scatter(0.2, mean(RT_fitslopes_left), '<k', 'filled')
% scatter(x*0.3, RT_fitslopes_leftLED, '<b')
% scatter(0.3, mean(RT_fitslopes_leftLED), '<b', 'filled')
% scatter(x*0.7, RT_fitslopes_right, '>k')
% scatter(0.7, mean(RT_fitslopes_right), '>k', 'filled')
% scatter(x*0.8, RT_fitslopes_rightLED, '>b')
% scatter(0.8, mean(RT_fitslopes_rightLED), '>b', 'filled')
% xlim([0 1])
% title('RT - fit slopes')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('RT')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('RT_fitslopes', '-dpdf') 
% 
% figure;
% hold on
% x = ones(1, length(RT_intercepts_left));
% scatter(x*0.2, RT_intercepts_left, '<k')
% scatter(0.2, mean(RT_intercepts_left), '<k', 'filled')
% scatter(x*0.3, RT_intercepts_leftLED, '<b')
% scatter(0.3, mean(RT_intercepts_leftLED), '<b', 'filled')
% scatter(x*0.7, RT_intercepts_right, '>k')
% scatter(0.7, mean(RT_intercepts_right), '>k', 'filled')
% scatter(x*0.8, RT_intercepts_rightLED, '>b')
% scatter(0.8, mean(RT_intercepts_rightLED), '>b', 'filled')
% xlim([0 1])
% title('RT - fit intercepts')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('RT')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('RT_intercepts', '-dpdf') 
% 
% %
% 
% figure;
% hold on
% x = ones(1, length(slopes_left));
% scatter(x*0.2, slopes_left, '<k')
% scatter(0.2, nanmean(slopes_left), '<k', 'filled')
% scatter(x*0.3, slopes_leftLED, '<b')
% scatter(0.3, nanmean(slopes_leftLED), '<b', 'filled')
% scatter(x*0.7, slopes_right, '>k')
% scatter(0.7, nanmean(slopes_right), '>k', 'filled')
% scatter(x*0.8, slopes_rightLED, '>b')
% scatter(0.8, nanmean(slopes_rightLED), '>b', 'filled')
% xlim([0 1])
% title('mean slopes')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('slope')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('mean_slopes', '-dpdf') 
% 
% figure;
% hold on
% x = ones(1, length(RTs_left));
% scatter(x*0.2, RTs_left, '<k')
% scatter(0.2, nanmean(RTs_left), '<k', 'filled')
% scatter(x*0.3, RTs_leftLED, '<b')
% scatter(0.3, nanmean(RTs_leftLED), '<b', 'filled')
% scatter(x*0.7, RTs_right, '>k')
% scatter(0.7, nanmean(RTs_right), '>k', 'filled')
% scatter(x*0.8, RTs_rightLED, '>b')
% scatter(0.8, nanmean(RTs_rightLED), '>b', 'filled')
% xlim([0 1])
% title('mean RTs')
% ax = gca;
% ax.XTick = [0.2, 0.3, 0.7, 0.8];
% ax.XTickLabels = {'L-Cntrl', 'L-LED', 'R-Cntrl', 'R-LED'};
% ylabel('RT')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('mean_RTs', '-dpdf')



%% fit stats

% [h_RT_fitslopes_left, p_RT_fitslopes_left] = ttest(RT_fitslopes_left);
% [h_RT_fitslopes_leftLED, p_RT_fitslopes_leftLED] = ttest(RT_fitslopes_leftLED);
% [h_RT_fitslopes_right, p_RT_fitslopes_right] = ttest(RT_fitslopes_right);
% [h_RT_fitslopes_rightLED, p_RT_fitslopes_rightLED] = ttest(RT_fitslopes_rightLED);
% 
% [h_RT_intercepts_left, p_RT_intercepts_left] = ttest(RT_intercepts_left);
% [h_RT_intercepts_leftLED, p_RT_intercepts_leftLED] = ttest(RT_intercepts_leftLED);
% [h_RT_intercepts_right, p_RT_intercepts_right] = ttest(RT_intercepts_right);
% [h_RT_intercepts_rightLED, p_RT_intercepts_rightLED] = ttest(RT_intercepts_rightLED);
% 
% names = ([{'p_intercepts'}; {'p_slopes'}; {'h_intercepts'}; {'h_slopes'}])
% left = [p_RT_intercepts_left; p_RT_fitslopes_left; h_RT_intercepts_left ;h_RT_fitslopes_left]
% leftLED = [p_RT_intercepts_leftLED; p_RT_fitslopes_leftLED; h_RT_intercepts_leftLED  ;h_RT_fitslopes_leftLED]
% right = [p_RT_intercepts_right; p_RT_fitslopes_right; h_RT_intercepts_right  ;h_RT_fitslopes_right]
% rightLED = [p_RT_intercepts_rightLED; p_RT_fitslopes_rightLED; h_RT_intercepts_rightLED ;h_RT_fitslopes_rightLED]
% 
% Ttest_results = table(names, left,leftLED, right, rightLED);
% 
% 
% %L and R ;, columns = side,  rows = stim mice are replicates
% RT_fitslopes = [[RT_fitslopes_left(:); RT_fitslopes_leftLED(:)] [RT_fitslopes_right(:) ;RT_fitslopes_rightLED(:)]];
% nmice = 5;
% [p_RT_fitslopes, tbl_RT_fitslopes, stats_RT_fitslopes] = anova2(RT_fitslopes, nmice);
% slope_p_stim = p_RT_fitslopes(1);
% slope_p_side = p_RT_fitslopes(2);
% slope_p_interaction = p_RT_fitslopes(3);
% 
% RT_intercepts = [[RT_intercepts_left(:); RT_intercepts_leftLED(:)] [RT_intercepts_right(:); RT_intercepts_rightLED(:)]];
% [p_RT_intercepts, tbl_RT_intercepts, stats_RT_intercepts] = anova2(RT_intercepts, nmice);
% intercepts_p_stim = p_RT_intercepts(1);
% intercepts_p_side = p_RT_intercepts(2);
% intercepts_p_interaction = p_RT_intercepts(3);
% 
% Factor = ([{'Side'}; {'Stim'}; {'Interaction'}]);
% slopes =  p_RT_fitslopes(:);
% intercepts =  p_RT_intercepts(:);
% table(Factor, slopes, intercepts)
% 
% RT.meanslopeLeft = [mean(RT_fitslopes_left), mean(RT_fitslopes_leftLED)];
% RT.meanslopeRight = [mean(RT_fitslopes_right), mean(RT_fitslopes_rightLED)];
% RT.meaninterceptLeft = [mean(RT_intercepts_left), mean(RT_intercepts_leftLED)];
% RT.meaninterceptRight = [mean(RT_intercepts_right), mean(RT_intercepts_rightLED)];



%%
% 
% dslope_fitslopes_left = slope_fitslopes_leftLED - slope_fitslopes_left;
% dslope_fitslopes_right = slope_fitslopes_rightLED - slope_fitslopes_right;
% 
% dslope_intercepts_left = slope_intercepts_leftLED - slope_intercepts_left;
% dslope_intercepts_right = slope_intercepts_rightLED - slope_intercepts_right;
% 
% dRT_fitslopes_left = RT_fitslopes_leftLED - RT_fitslopes_left;
% dRT_fitslopes_right = RT_fitslopes_rightLED - RT_fitslopes_right;
% 
% dRT_intercepts_left = RT_intercepts_leftLED - RT_intercepts_left;
% dRT_intercepts_right = RT_intercepts_rightLED - RT_intercepts_right;
% 
% 
% figure;
% hold on
% x = categorical({'Left', 'Right'});
% x = reordercats(x,{'Left', 'Right'});
% y = [dslope_fitslopes_left; dslope_fitslopes_right ];
% b = bar(x, y);
% ylabel('delta slope')
% title('Slope - Change in fit slopes')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('Slope_delta_fitslope', '-dpdf')
% 
% figure;
% hold on
% x = categorical({'Left', 'Right'});
% x = reordercats(x,{'Left', 'Right'});
% y = [dslope_intercepts_left; dslope_intercepts_right];
% b = bar(x, y);
% ylabel('delta intercept')
% title('Slope - Change in fit intercepts')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('Slope_delta_intercept', '-dpdf')
% 
% figure;
% hold on
% x = categorical({'Left', 'Right'});
% x = reordercats(x,{'Left', 'Right'});
% y = [dRT_fitslopes_left; dRT_fitslopes_right ];
% b = bar(x, y);
% ylabel('delta slope')
% title('RT - Change in fit slopes')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('RT_delta_fitslope', '-dpdf')
% 
% figure;
% hold on
% x = categorical({'Left', 'Right'});
% x = reordercats(x,{'Left', 'Right'});
% y = [dRT_intercepts_left; dRT_intercepts_right];
% b = bar(x, y);
% ylabel('delta intercept')
% title('RT - Change in fit intercepts')
% 
% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\figures\test\')
% print('RT_delta_intercept', '-dpdf')




