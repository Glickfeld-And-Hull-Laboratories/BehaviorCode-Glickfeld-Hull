%% data concatenation
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\';
load(fullfile(data_path, 'CDMouseList.mat'))

% mouse_names_PV_contrast = ['i419'; 'i548'; 'i565'; 'i578'; 'i593']
% temp = struct();

% All Curve stats
% for m =1:length(mouse_names_PV_contrast)
%     mouse = allmice(m).name
%     load(fullfile(data_path,[mouse 'goodData.mat']));
%     temp(m) = input;
% end

load(fullfile(data_path,['i565' 'goodData.mat']));

% input = concatenateDataBlocks(temp); % Using input from "mouse"goodData.mat - CLM
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
        stimtime = double(input.stimTimestampMs{trN});
        stimTime= stimtime+delay;
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
            return
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
                                  rc_empty_control = rc_empty_control+1
                              else
                                  rc_empty_LED = rc_empty_LED+1
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  rc_full_control = rc_full_control+1
                              else
                                  rc_full_LED = rc_full_LED+1
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end  
                      elseif iloc==1 && ibehavior==1 %right incorrect (pos values)
                          ifirst_qval = find(iqvals >=ldt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  ri_empty_control = ri_empty_control+1
                              else
                                  ri_empty_LED = ri_empty_LED+1
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  ri_full_control = ri_full_control+1
                              else
                                  ri_full_LED = ri_full_LED+1
                              end
                              temptime_itrial = ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end                      
                      elseif iloc==2 && ibehavior==2 %left correct (pos values)
                          ifirst_qval = find(iqvals >=ldt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  lc_empty_control = lc_empty_control+1
                              else
                                  lc_empty_LED = lc_empty_LED+1
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  lc_full_control = lc_full_control+1
                              else
                                  lc_full_LED = lc_full_LED+1
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end
                      elseif iloc==2 && ibehavior==1 %left incorrect (neg values)
                          ifirst_qval = find(iqvals <=rdt,1,'first');
                          if isempty(ifirst_qval)
                              if icontrol == 1
                                  li_empty_control = li_empty_control+1
                              else
                                  li_empty_LED = li_empty_LED+1
                              end
%                               plot(iqvals)
                              continue
                          else
                              if icontrol == 1
                                  li_full_control = li_full_control+1
                              else
                                  li_full_LED = li_full_LED+1
                              end
                              temptime_itrial =ifirst_qval;
                              temptime = [temptime temptime_itrial];
                          end   
                      end
                      
                   rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))=temptime; %times (where quad reached threshold) across trials matching conditions
                   rxntime_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(temptime); % average times (where quad reached threshold) matching conditions
                   
                   rxnslopewindowstart=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))-50;
                   rxnslopewindowend=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))+400; % Why + 400 and not +50?
                   
                   rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime))=qvals_temp_mean(temptime+8000);    
                   rxnqval_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime)));

                   qvals_adj_window = qvals_adj(8000:18000, :);
                   rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime))=abs((qvals_adj(rxnslopewindowend+8000)-qvals_adj(rxnslopewindowstart+8000)))/450;
                   rxnslope_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   
%                    rxnslope_succ_fail_avg(icontrol,:,iloc,irat,1)= mean(nanmean(rxnslope(icontrol,:,iloc,irat,1:length(temptime)))); % combined correct, incorrect trials
                   
                   rxntime_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxnslope_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));                 
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

percent_rc_control = rc_empty_control / (rc_empty_control + rc_full_control)
percent_ri_control = ri_empty_control/ (ri_empty_control + ri_full_control)
percent_lc_control = lc_empty_control/ (lc_empty_control + lc_full_control)
percent_li_control = li_empty_control/ (li_empty_control + li_full_control)

percent_rc_LED = rc_empty_LED / (rc_empty_LED + rc_full_LED)
percent_ri_LED = ri_empty_LED/ (ri_empty_LED + ri_full_LED)
percent_lc_LED = lc_empty_LED/ (lc_empty_LED + lc_full_LED)
percent_li_LED = li_empty_LED/ (li_empty_LED + li_full_LED)

percent.rc_control = percent_rc_control
percent.ri_control = percent_ri_control
percent.lc_control = percent_lc_control
percent.li_control= percent_li_control
percent.rc_LED= percent_rc_LED
percent.ri_LED= percent_ri_LED
percent.lc_LED= percent_lc_LED
percent.li_LED= percent_li_LED
%% Find %correct by condition

%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1 --- 4 is highest contrast

for icontrol=1:nledcond
   for iloc = 1:nstimlocation
       for irat = 1:nratios
%            correct = sum(~isnan(sorted_trials(icontrol,2,iloc,irat,:)));
%            incorrect = sum(~isnan(sorted_trials(icontrol,1,iloc,irat,:)));
           correct = length(sorted_trials{icontrol,2,iloc,irat});
           incorrect = length(sorted_trials{icontrol,1,iloc,irat});
           total = correct + incorrect;
           prct_corr(icontrol, iloc, irat) = correct / total;
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
    
    %error
    
    mscC_left_error(i) = rxnslope_error(1,2,2,i);
    mscC_right_error(i) = rxnslope_error(1,2,1,i);
    
    mscLED_left_error(i) = rxnslope_error(2,2,2,i);
    mscLED_right_error(i) = rxnslope_error(2,2,1,i);
    
    msiC_left_error(i) = rxnslope_error(1,1,2,i);
    msiC_right_error(i) = rxnslope_error(1,1,1,i);
    
    msiLED_left_error(i) = rxnslope_error(2,1,2,i);
    msiLED_right_error(i) = rxnslope_error(2,1,1,i);
    
    %time to thresh
    
    mtcC_left(i) = rxntime_avg(1,2,2,i);
    mtcC_right(i) = rxntime_avg(1,2,1,i);
    
    mtcLED_left(i) = rxntime_avg(2,2,2,i);
    mtcLED_right(i) = rxntime_avg(2,2,1,i);
    
    mtiC_left(i) = rxntime_avg(1,1,2,i);
    mtiC_right(i) = rxntime_avg(1,1,1,i);
    
    mtiLED_left(i) = rxntime_avg(2,1,2,i);
    mtiLED_right(i) = rxntime_avg(2,1,1,i);
    
    %combined correct/incorrect slope
    
    meanslope_all__control_left(i) = rxnslope_succ_fail_avg(1,1,2,i);
    meanslope_all__control_right(i) = rxnslope_succ_fail_avg(1,1,1,i);
    
    meanslope_all__LED_left(i) = rxnslope_succ_fail_avg(2,1,2,i);
    meanslope_all__LED_right(i) = rxnslope_succ_fail_avg(2,1,1,i);
    
    %prct correct
    
    prct_corr_left_control(i) = prct_corr(1,2,i);
    prct_corr_right_control(i) = prct_corr(1,1,i);
    prct_corr_left_LED(i) = prct_corr(2,2,i);
    prct_corr_right_LED(i) = prct_corr(2,1,i);
end

%slope
meanslope_correct_control=horzcat(fliplr(mscC_left),mscC_right);
meanslope_correct_LED=horzcat(fliplr(mscLED_left),mscLED_right);

meanslope_incorrect_control=horzcat(fliplr(msiC_left),msiC_right);
meanslope_incorrect_LED=horzcat(fliplr(msiLED_left),msiLED_right);

meanslope_all_control = horzcat(fliplr(meanslope_all__control_left),meanslope_all__control_right);
meanslope_all_LED = horzcat(fliplr(meanslope_all__LED_left),meanslope_all__LED_right);

%error

meanslope_correct_control_error=horzcat(fliplr(mscC_left_error),mscC_right_error);
meanslope_correct_LED_error=horzcat(fliplr(mscLED_left_error),mscLED_right_error);

meanslope_incorrect_control_error=horzcat(fliplr(msiC_left_error),msiC_right_error);
meanslope_incorrect_LED_error=horzcat(fliplr(msiLED_left_error),msiLED_right_error);



%time

meantime_correct_control=horzcat(fliplr(mtcC_left),mtcC_right);
meantime_correct_LED=horzcat(fliplr(mtcLED_left),mtcLED_right);

meantime_incorrect_control=horzcat(fliplr(mtiC_left),mtiC_right);
meantime_incorrect_LED=horzcat(fliplr(mtiLED_left),mtiLED_right);

xright=[c4 c3 c2 c1]; xleft=[1/c1 1/c2 1/c3 1/c4];
x=horzcat(xleft,xright);
xstr=sprintf('%s ratio (R/L)',taskstr);



figure;
scatter(x ,meanslope_correct_control, 'k')
hold on
scatter(x ,meanslope_correct_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
set(gca,'box','off')
xlabel(xstr), ylabel('Slope'); 


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

               
               
               
               