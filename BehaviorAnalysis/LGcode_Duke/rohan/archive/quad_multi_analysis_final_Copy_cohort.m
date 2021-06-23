%% data concatenation
data_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\';
load(fullfile(data_path, 'CDMouseList.mat'))
%%
mouse_names_PV_contrast = ['i419'; 'i548'; 'i565'; 'i578'; 'i593']


% All Curve stats
for m =1:length(mouse_names_PV_contrast)
    mouse = allmice(m).name
    load(fullfile(data_path,[mouse 'goodData.mat']));
    
    
   tempStruct(m).trialOutcomeCell = input.trialOutcomeCell;
   tempStruct(m).tDecisionTimeMs = input.tDecisionTimeMs;
   tempStruct(m).tBlock2TrialNumber = input.tBlock2TrialNumber;
   tempStruct(m).tGratingEccentricityDeg = input.tGratingEccentricityDeg;
   tempStruct(m).feedbackMotionSensitivity = input.feedbackMotionSensitivity;
   tempStruct(m).tGratingContrast = input.tGratingContrast;
   tempStruct(m).dGratingContrast  = input.dGratingContrast;
   tempStruct(m).tLeftTrial = input.tLeftTrial;
   tempStruct(m).tThisTrialStartTimeMs = input.tThisTrialStartTimeMs;
   tempStruct(m).stimTimestampMs = input.stimTimestampMs;
   tempStruct(m).tDecisionTimeMs = input.tDecisionTimeMs;
   tempStruct(m).block2TrialLaserPowerMw = input.block2TrialLaserPowerMw;
   tempStruct(m).quadratureTimesUs  = input.quadratureTimesUs;
   tempStruct(m).quadratureValues  = input.quadratureValues;
 
    
end

% load(fullfile(data_path,['i565' 'goodData.mat']));

input = concatenateDataBlocks(tempStruct); % Using input from "mouse"goodData.mat - CLM
%% 


% cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
% savedate=input.saveTime; 
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell))); %%? Why 18001 (total time window? 8000 + 10000)? Why use uint16 (changes class to unit16 - 16-bit intergers)? - CLM
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
SIx = double(strcmp(input.trialOutcomeCell,'success'));
FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
block2=celleqel2mat_padded(input.tBlock2TrialNumber);
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg); % 
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity(1); rdt=-1*ldt; % ldt/rdt = left/right decision threshold - CLM 
tContrast=celleqel2mat_padded(input.tGratingContrast); % comment out for size discrim
unqtargets=unique(tContrast);
dContrast=celleqel2mat_padded(input.dGratingContrast); % comment out for size discrim
%tContrast=celleqel2mat_padded(input.tGratingDiameterDeg); 
%dContrast=celleqel2mat_padded(input.dGratingDiameterDeg);
contrastratio=tContrast./dContrast;
contrastratio=round(contrastratio,3);
uncontrasts=unique(contrastratio);
c1=uncontrasts(4); c2=uncontrasts(3); c3=uncontrasts(2); c4=uncontrasts(1);
tLeftTrial=celleqel2mat_padded(input.tLeftTrial);
numtrials=length(tLeftTrial);
thistrialtime=celleqel2mat_padded(input.tThisTrialStartTimeMs);
stimtime=celleqel2mat_padded(input.stimTimestampMs);
dectime=celleqel2mat_padded(input.tDecisionTimeMs);
ledpwr=input.block2TrialLaserPowerMw;
for trN = 1:length(input.trialOutcomeCell)-1 
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000); %  Times of each registered wheel movement starting from current and ending of next trial. Microseconds to milliseconds
        qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]); % Wheel positions from current trial to end of next
        stimTime = double(input.stimTimestampMs{trN});
        qTimes_zero = qTimes-stimTime; % Wheel movement before (-) and after (+) stim appears
        qVals = qVals-qVals(1); % Start quad values from 0
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000); % Index of times where wheel movement falls in time window
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind); % Time of wheel movement relative, to stimulus appearance, within time window
            qVals_sub = qVals(time_ind); % Wheel position in time window
            qTimes_temp = qTimes(time_ind); % Time of raw (non-aligned) wheel movement in time window
            rep_ind = find(diff(qTimes_sub)==0); % ? Returns index where the difference between adjacent elements of q_times_sub is 0. Repeated time
            qTimes_sub(rep_ind) = []; % removes rep_ind (repeated index) from vector
            qVals_sub(rep_ind) = [];
            qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000; % Time window
            qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)'; % Not sure why this is being done
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)'; % interpolated movement
            if input.tDecisionTimeMs{trN} < 10000
                if isnan(qVals_final(8000,trN))
                    qVals_final(8000,trN) = qVals_final(find(~isnan(qVals_final(:,trN)),1,'first'),trN); % find first movment so can align to 0
                end
                qVal_off = qVals_final(:,trN)-qVals_final(8000,trN); %align final quad values to 0, currently not used
                qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN); % currently not used...
            end
        else
            return
        end; end; end
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
sorted_trials = nan(nledcond, ncorrectcond, nstimlocation, nratios, numtrials);
for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear tempcon tempbehavesuccess tempbehavefailure tempbehave temploc temprat temp
               tempcon = find(block2 == icontrol-1);
               if ibehavior == 1
                   tempbehave = find(FIx == 1);
               elseif ibehavior == 2
                   tempbehave = find(SIx == 1);
               end
               temploc = find(tLeftTrial == iloc-1);
               temprat = find(contrastratio == uncontrasts(irat));
               temp =intersect(intersect(intersect(tempcon, tempbehave), temploc), temprat); % indexes of trials where these conditions are aligned
               sorted_trials(icontrol, ibehavior, iloc, irat, 1:length(temp)) = temp; % sorted these aligned conditons into matrix by condition
           end; end; end; end
%%
bluecolors=brewermap(11, '*Blues');
redcolors=brewermap(11, '*Reds');
rxntime=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxntime_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxntime_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope_error=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear temp qvals_adj qvals_temp_mean time_adj temptime_itrial temptime
                if icontrol==1
                       usecolor=bluecolors;
                   elseif icontrol==2
                       usecolor=redcolors; 
                    end
               temp = squeeze(sorted_trials(icontrol,ibehavior,iloc,irat,:));               
               temp=temp(find(~isnan(temp)));
               qvals_adj=qVals_final(:,temp)-qVals_final(8000,temp);               
               qvals_temp_mean=nanmean(qvals_adj,2);  % mean qVals of trials that share same conditions
               temptime = [ ];
                   for itrial=1:length(temp)
                      if iloc==1 && ibehavior==2 %right correct
                          if isempty(min(find(qvals_adj(8000:18000,itrial)<=rdt)))
                              continue
                          else
                            temptime_itrial =min(find(qvals_adj(8000:18000,itrial)<=rdt)); % removed (itrial) after temptime1 - CLM
                            temptime = [temptime temptime_itrial];
                          end  
                      elseif iloc==1 && ibehavior==1 %right incorrect
                          if isempty(min(find(qvals_adj(8000:18000,itrial)>=ldt)))
                              continue
                          else
                          temptime_itrial =min(find(qvals_adj(8000:18000,itrial)>=ldt));
                          temptime = [temptime temptime_itrial];
                          end                      
                      elseif iloc==2 && ibehavior==2 %left correct
                          if isempty(min(find(qvals_adj(8000:18000,itrial)>=ldt)))
                              continue
                          else
                          temptime_itrial =min(find(qvals_adj(8000:18000,itrial)>=ldt));
                          temptime = [temptime temptime_itrial];
                          end
                      elseif iloc==2 && ibehavior==1 %left incorrect
                          if isempty(min(find(qvals_adj(8000:18000,itrial)<=rdt)))
                              continue
                          else
                          temptime_itrial =min(find(qvals_adj(8000:18000,itrial)<=rdt));
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
                   rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime))=abs((qvals_adj_window(rxnslopewindowend)-qvals_adj_window(rxnslopewindowstart)))/450;
                   rxnslope_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   
                   rxnslope_succ_fail_avg(icontrol,:,iloc,irat,1)= mean(nanmean(rxnslope(icontrol,:,iloc,irat,1:length(temptime)))); % combined correct, incorrect trials
                   
                   rxntime_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxnslope_error(icontrol,ibehavior,iloc,irat,1)=nanstd(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));                 
                   end
                
               if iloc==1 && ibehavior==2 && icontrol==1
                   subplot(2,4,1)
                   hold on
                   title('Right Correct Control Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1); %
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==2 && icontrol==2
                   subplot(2,4,2)
                   hold on
                   title('Right Correct LED Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==2 && icontrol==1
                   subplot(2,4,3)
                   hold on
                   title('Left Correct Control Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                    line([0 0], [-500 500],'Color','k')
                   line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==2 && icontrol==2
                   subplot(2,4,4)
                   hold on
                   title('Left Correct LED Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==1 && icontrol==1
                   subplot(2,4,5)
                   hold on
                   title('Right Incorrect Control Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==1 && icontrol==2
                   subplot(2,4,6)
                   hold on
                   title('Right Incorrect LED Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==1 && icontrol==1
                   subplot(2,4,7)
                   hold on
                   title('Left Incorrect Control Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k')
                   line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==1 && icontrol==2
                   subplot(2,4,8)
                   hold on
                   title('Left Incorrect LED Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   checkthresh = rxntime_avg(icontrol,ibehavior,iloc,irat,1);
                   if isnan(checkthresh) %% Go back and check left incorrect values. May not want to use this method.
                       continue
                   else
                       plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   end
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
                   end;
               hold on
           end; end; end; end               
%% 
% scatterplot %correct (left/right) vs. slope (by contrast ratio) (All trials, Control vs LED
% same plot

% scatterrplot %incorrect (left/right) vs. slope (by contrast ratio) (Incorrect trials, Control vs LED
% same plot

% if input.doSizeDiscrim==1
%     taskstr='size';
% if input.doContrastDiscrim==1
    taskstr='Contrast';
% end
numignores=length(Iix);
num_includedtrials=numtrials-numignores;
numcorrect=sum(SIx); numincorrect=sum(FIx);
numcontrols=numtrials-sum(block2); numled=sum(block2);



 % load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_discrim\data\final\i565_pctrightInfo.mat')
 



% meanRT_correct=horzcat(fliplr(meanRT_correct_left), meanRT_correct_right); % mean RT for all correct (Control vs LED) 
% meanRT_incorrect=horzcat(fliplr(meanRT_incorrect_left), meanRT_incorrect_right); % mean RT for all incorrect
% SEMRT_correct=horzcat(fliplr(SEMRT_correct_left), SEMRT_correct_right);
% SEMRT_incorrect=horzcat(fliplr(SEMRT_incorrect_left), SEMRT_incorrect_right);

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

% 
% meanslope_incorrect=horzcat(fliplr(meanslope_incorrect_left_control), meanslope_incorrect_right_control);
% SEMslope_correct=horzcat(fliplr(SEMslope_correct_left_control),SEMslope_correct_right_control);
% SEMslope_incorrect=horzcat(fliplr(SEMslope_incorrect_left_control),SEMslope_incorrect_right_control);

xright=[c4 c3 c2 c1]; xleft=[1/c1 1/c2 1/c3 1/c4];
x=horzcat(xleft,xright);
xstr=sprintf('%s ratio (R/L)',taskstr);

figure(2)
scatter(x, meanslope_correct_control, 'k')
% errorbar(x, meanslope_correct_control, meanslope_correct_control_error,'ok')
hold on
scatter(x, meanslope_correct_LED, 'b')
% errorbar(x, meanslope_correct_LED, meanslope_correct_LED_error,'ob')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
xlabel(xstr); ylabel('Quadrature trace slope'); 
titlestr=sprintf('Avg quad response slope by contrast ratio - Correct');
title(titlestr)
% legend('Control', 'LED')


figure(3)
scatter(x, meanslope_incorrect_control, 'k')
hold on
scatter(x, meanslope_incorrect_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
xlabel(xstr); ylabel('Quadrature trace slope'); 
titlestr=sprintf('Avg quad response slope by contrast ratio - Incorrect');
title(titlestr)

figure(4)
scatter(x, meantime_correct_control, 'k')
hold on
scatter(x, meantime_correct_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
xlabel(xstr); ylabel('Response Time'); 
titlestr=sprintf('Avg response time by contrast ratio - Correct');
title(titlestr)
% legend('Control', 'LED')


figure(5)
scatter(x, meantime_incorrect_control, 'k')
hold on
scatter(x, meantime_incorrect_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
xlabel(xstr); ylabel('Response Time'); 
titlestr=sprintf('Avg response time by contrast ratio - Incorrect');
title(titlestr)

figure(6)
scatter3(x,pct_Corr_all(1,:),meanslope_all_control, 'k')
hold on
scatter3(x,pct_Corr_all(2,:),meanslope_all_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
xlabel(xstr); ylabel('% Correct'), zlabel('Slope'); 
title('Avg quad slope by %Correct and Contrast Ratio');


figure(7)
scatter(x ,meanslope_all_control, 'k')
hold on
scatter(x ,meanslope_all_LED, 'b')
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlim([0 102]);
set(gca,'box','off')
xlabel(xstr), ylabel('Slope'); 

% subaxis(2,2,1,'SpacingVert',0.02,'MR',0.01); 
% hold on
% pc_control=plot(x,meanRT_correct(1,:),'-blueo','MarkerFaceColor','blue')
% pi_control=plot(x,meanRT_incorrect(1,:),'-redo','MarkerFaceColor','red')
% errorbar(x,meanRT_correct(1,:), SEMRT_correct(1,:), 'blueo')
% errorbar(x,meanRT_incorrect(1,:), SEMRT_incorrect(1,:), 'redo')
% xlim([0 102]); ylim([0 8000])
% ylabel('Milliseconds after stim on'); 
% ax = gca; ax.XScale = 'log'
% set(gca,'XTickLabel', [])
% set(gca,'box','off')
% legend('Control correct', 'Control incorrect')
% titlestr=sprintf('Average subjective reaction time by R/L %s ratio' , taskstr);
% title(titlestr);

% subaxis(2,2,3,'SpacingVert',0.02,'MR',0.01); 
% hold on
% pc_led=plot(x,meanRT_correct(2,:),'-cyano','MarkerFaceColor','cyan')
% pi_led=plot(x,meanRT_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
% errorbar(x,meanRT_correct(2,:), SEMRT_correct(2,:), 'cyano')
% errorbar(x,meanRT_incorrect(2,:), SEMRT_incorrect(2,:), 'magentao')
% xlabel(xstr); ylabel('Milliseconds after stim on'); 
% xlim([0 102]);ylim([0 8000])
% ax = gca; ax.XScale = 'log'
% set(gca,'box','off')
% legend('LED correct', 'LED incorrect')
% 
% subaxis(2,2,2,'SpacingVert',0.02,'MR',0.01); 
% hold on
% pcc_acttime=plot(x,meanacttime_correct(1,:),'-blueo','MarkerFaceColor','blue')
% pci_acttime=plot(x,meanacttime_incorrect(1,:),'-redo','MarkerFaceColor','red')
% errorbar(x,meanacttime_correct(1,:),SEMacttime_correct(1,:),'blueo','LineStyle','none');
% errorbar(x,meanacttime_incorrect(1,:), SEMacttime_incorrect(1,:),'redo','LineStyle','none');
% xlim([0 102]);ylim([0 8000])
% ax = gca; ax.XScale = 'log'
% set(gca,'XTickLabel',[],'YTickLabel', [])
% set(gca,'box','off')
% titlestr=sprintf('Average true decision time by R/L %s ratio' , taskstr);
% title(titlestr);
% 
% subaxis(2,2,4,'SpacingVert',0.02,'MR',0.01); 
% hold on
% pledi_acttime=plot(x,meanacttime_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
% pledc_acttime=plot(x,meanacttime_correct(2,:),'-cyano','MarkerFaceColor','cyan')
% errorbar(x,meanacttime_correct(2,:),SEMacttime_correct(2,:),'cyano','LineStyle','none');
% errorbar(x,meanacttime_incorrect(2,:), SEMacttime_incorrect(2,:),'magentao','LineStyle','none');
% xlim([0 102]);ylim([0 8000])
% xlabel(xstr);  
% set(gca,'YTickLabel', [])
% ax = gca; ax.XScale = 'log'
% set(gca,'box','off')
% 
% figure(3)
% subaxis(1,2,1,'SpacingHor',0.02,'MR',0.01); 
% hold on
% xlim([0 102]); ylim([0 0.2])
% ax = gca
% ax.XScale = 'log'
% set(gca,'box','off')
% xlabel(xstr); ylabel('Absolute value of quadrature trace slope'); 
% titlestr=sprintf('Average quadrature response trace slope by R/L %s ratio.', taskstr);
% title(titlestr)
% pcc_rxnslope=plot(x,meanslope_correct(1,:),'-blueo','MarkerFaceColor','blue')
% pci_rxnslope=plot(x,meanslope_incorrect(1,:),'-redo','MarkerFaceColor','red')
% errorbar(x,meanslope_correct(1,:),SEMslope_correct(1,:),'blueo','LineStyle','none');
% errorbar(x,meanslope_incorrect(1,:), SEMslope_incorrect(1,:),'redo','LineStyle','none');
% legend('Control Correct', 'Control Incorrect')
% 
% subaxis(1,2,2,'SpacingHor',0.02,'MR',0.01); 
% hold on
% xlim([0 102]); ylim([0 0.2])
% ax=gca
% ax.XScale = 'log'
% xlabel(xstr);
% set(gca, 'YTickLabel',[])
% set(gca,'box','off')
% titlestr=sprintf('Average quadrature response trace slope by R/L %s ratio.' , taskstr);
% title(titlestr)
% pledc_rxnslope=plot(x,meanslope_correct(2,:),'-cyano','MarkerFaceColor','cyan')
% pledi_rxnslope=plot(x,meanslope_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
% errorbar(x,meanslope_correct(2,:),SEMslope_correct(2,:),'cyano','LineStyle','none');
% errorbar(x,meanslope_incorrect(2,:), SEMslope_incorrect(2,:),'magentao','LineStyle','none');
% legend('LED Correct', 'LED Incorrect')

%%
clearvars -except rxntime_correct_right rxntime_correct_left rxntime_incorrect_right rxntime_incorrect_left...
    rxnslope_correct_right rxnslope_correct_left rxnslope_incorrect_right rxnslope_incorrect_left...
    meanslope_correct_right meanslope_correct_left meanslope_incorrect_right meanslope_incorrect_left...
    meanRT_correct_right meanRT_correct_left meanRT_incorrect_right meanRT_incorrect_left...
    SEMslope_correct_right SEMslope_correct_left SEMslope_incorrect_right SEMslope_incorrect_left...
    SEMRT_correct_right SEMRT_correct_left SEMRT_incorrect_right SEMRT_incorrect_left ...
    acttime_correct_right acttime_correct_left acttime_incorrect_right...
    acttime_incorrect_left SEMacttime_correct_right SEMacttime_correct_left SEMacttime_incorrect_right...
    SEMacttime_incorrect_left meanacttime_correct_left meanacttime_correct_right meanacttime_incorrect_right...
    meanacttime_incorrect_left

               
               
               
               