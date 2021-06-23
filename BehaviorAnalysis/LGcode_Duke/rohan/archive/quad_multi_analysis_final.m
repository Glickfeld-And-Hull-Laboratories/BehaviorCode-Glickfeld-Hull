%% data concatenation
% input = concatenateDataBlocks(temp); % Using input from "mouse"goodData.mat - CLM
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
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity; rdt=-1*ldt; % ldt/rdt = left/right decision threshold - CLM 
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
for trN = 1:length(input.trialOutcomeCell)-1 %not clear on why -1 - CLM
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000); %  Times of each registered wheel movement on each trial?
        qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]); % Wheel positions from current trial to end of next
        stimTime = double(input.stimTimestampMs{trN});
        qTimes_zero = qTimes-stimTime; % Wheel movement before (-) and after (+) stim appears
        qVals = qVals-qVals(1); % ???
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000); % Index where wheel movement falls in time window?
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind); % Time of wheel movement relative, to stimulus appearance, within time window
            qVals_sub = qVals(time_ind); % Wheel position in time window
            qTimes_temp = qTimes(time_ind); % Time of wheel movement in time window
            rep_ind = find(diff(qTimes_sub)==0); % ? Returns index where the difference between adjacent elements of q_times_sub is 0...why? Repeated times?
            qTimes_sub(rep_ind) = []; % removes rep_ind
            qVals_sub(rep_ind) = [];
            qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000; % Time window
            qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)'; % ?? 1-D data interpolation?? (Value approximation/ curve fitting)
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)'; % ??
            if input.tDecisionTimeMs{trN} < 10000
                if isnan(qVals_final(8000,trN))
                    qVals_final(8000,trN) = qVals_final(find(~isnan(qVals_final(:,trN)),1,'first'),trN);
                end
                qVal_off = qVals_final(:,trN)-qVals_final(8000,trN);
                qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN);
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
nstimlocation=length(unique(eccentricity));
nratios=length(uncontrasts);
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1
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
               temp =intersect(intersect(intersect(tempcon, tempbehave), temploc), temprat);
               sorted_trials(icontrol, ibehavior, iloc, irat, 1:length(temp)) = temp;
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
               clear temp qvals_adj qvals_temp_mean time_adj temptime temptime1
                if icontrol==1
                       usecolor=bluecolors;
                   elseif icontrol==2
                       usecolor=redcolors; 
                    end
               temp = squeeze(sorted_trials(icontrol,ibehavior,iloc,irat,:));               
               temp=temp(find(~isnan(temp)));
               qvals_adj=qVals_final(:,temp)-qVals_final(8000,temp);               
               qvals_temp_mean=nanmean(qvals_adj,2);
                   for itrial=1:length(temp)
                      if iloc==1 && ibehavior==2 %right correct
                          temptime1(itrial) =min(find(qvals_adj(8000:18000,itrial)<=rdt)); 
                          if isempty(temptime1)
                              temptime(itrial)=min(find(qvals_adj(8000:18000)>=ldt));
                          else
                              temptime(itrial)=temptime1; end      
                      elseif iloc==1 && ibehavior==1 %right incorrect
                          temptime1(itrial) =min(find(qvals_adj(8000:18000,itrial)>=ldt)); 
                          if isempty(temptime1)
                              temptime(itrial)=min(find(qvals_adj(8000:18000)<=rdt));
                          else
                              temptime(itrial)=temptime1; end                            
                      elseif iloc==2 && ibehavior==2 %left correct
                          temptime1(itrial) =min(find(qvals_adj(8000:18000,itrial)>=ldt)); 
                          if isempty(temptime1)
                              temptime(itrial)=min(find(qvals_adj(8000:18000)<=rdt));
                          else 
                              temptime(itrial)=temptime1; end       
                      elseif iloc==2 && ibehavior==1 %left incorrect
                          temptime1(itrial) =min(find(qvals_adj(8000:18000,itrial)<=rdt));
                          if isempty(temptime1)
                              temptime(itrial)=min(find(qvals_adj(8000:18000)>=ldt));
                          else
                              temptime(itrial)=temptime1; end      
                      end                      
                   rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))=temptime;
                   rxntime_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(temptime);
                   rxnslopewindowstart=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))-50;
                   rxnslopewindowend=rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))+400;
                   rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime))=qvals_temp_mean(temptime+8000);    
                   rxnqval_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime))=abs((qvals_adj(rxnslopewindowend)-qvals_adj(rxnslopewindowstart)))/450;
                   rxnslope_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxntime_error=nanstd(rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime)));
                   rxnslope_error=nanstd(rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime)));                 
                  end
               if iloc==1 && ibehavior==2 && icontrol==1
                   subplot(2,4,1)
                   title('Right Correct Control Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==2 && icontrol==2
                   subplot(2,4,2)
                   title('Right Correct LED Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==2 && icontrol==1
                   subplot(2,4,3)
                   title('Left Correct Control Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                    line([0 0], [-500 500],'Color','k')
                   line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==2 && ibehavior==2 && icontrol==2
                   subplot(2,4,4)
                   title('Left Correct LED Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==1 && icontrol==1
                   subplot(2,4,5)
                   title('Right Incorrect Control Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               elseif iloc==1 && ibehavior==1 && icontrol==2
                   subplot(2,4,6)
                   title('Right Incorrect LED Trials')
                   ylim([-100 200]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
                   elseif iloc==2 && ibehavior==1 && icontrol==1
                   subplot(2,4,7)
                   title('Left Incorrect Control Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                    line([0 0], [-500 500],'Color','k')
                   line([800 800], [-500 500],'Color','k', 'LineStyle', '--')                   
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
                   elseif iloc==2 && ibehavior==1 && icontrol==2
                   subplot(2,4,8)
                   title('Left Incorrect LED Trials')
                   ylim([-200 100]); xlim([-500 6000]);
                   plot(qTimes_final, qvals_temp_mean,'Color',usecolor(irat*2,:));
                   line([0 0], [-500 500],'Color','k'); line([800 800], [-500 500],'Color','k', 'LineStyle', '--')
                   scatter(rxntime_avg(icontrol,ibehavior,iloc,irat,1), rxnqval_avg(icontrol,ibehavior,iloc,irat,1),'x','MarkerEdgeColor', usecolor(irat*2,:), 'LineWidth', 1.5);
               end
               hold on
           end; end; end; end               
%%
% if input.doSizeDiscrim==1
%     taskstr='size';
if input.doContrastDiscrim==1
    taskstr='contrast';
end
numignores=length(Iix);
num_includedtrials=numtrials-numignores;
numcorrect=sum(SIx); numincorrect=sum(FIx);
numcontrols=numtrials-sum(block2); numled=sum(block2);
for icontrol=1:nledcond  
    for iloc=1:nstimlocation    
        for irat=1:nratios            
            clear temptime_correct temptime_incorrect tempslope_correct tempslope_incorrect                
           temptime_correct=rxntime(icontrol,2,iloc,irat);
           temptime_incorrect=rxntime(icontrol,1,iloc,irat); 
           tempslope_correct=rxnslope(icontrol,2,iloc,irat);
           tempslope_incorrect=rxnslope(icontrol,1,iloc,irat);
           if ~isnan(squeeze(sorted_trials(icontrol,2,iloc,irat)))            
               clear acttime_correct
               acttime_correct=dectime(squeeze(sorted_trials(icontrol,2,iloc,irat)));
           end
           if ~isnan(squeeze(sorted_trials(icontrol,1,iloc,irat)));  
               clear acttime_incorrect
               acttime_incorrect=dectime(squeeze(sorted_trials(icontrol,1,iloc,irat))); 
           end
%            if iloc==1               
%                 rxntime_correct_right(icontrol,irat,ndatarun)=temptime_correct;
%                 rxntime_incorrect_right(icontrol,irat,ndatarun)=temptime_incorrect;  
%                 rxnslope_correct_right(icontrol,irat,ndatarun)=tempslope_correct;
%                 rxnslope_incorrect_right(icontrol,irat,ndatarun)=tempslope_incorrect;
%                 acttime_correct_right(icontrol,irat,ndatarun)=acttime_correct;
%                 acttime_incorrect_right(icontrol,irat,ndatarun)=acttime_incorrect;
%            elseif iloc==2
%                 rxntime_correct_left(icontrol,irat,ndatarun)=temptime_correct;
%                 rxntime_incorrect_left(icontrol,irat,ndatarun)=temptime_incorrect;  
%                 rxnslope_correct_left(icontrol,irat,ndatarun)=tempslope_correct;
%                 rxnslope_incorrect_left(icontrol,irat,ndatarun)=tempslope_incorrect;
%                 acttime_correct_left(icontrol,irat,ndatarun)=acttime_correct;
%                 acttime_incorrect_left(icontrol,irat,ndatarun)=acttime_incorrect;
%            end
        end; end; end;

meanRT_correct=horzcat(fliplr(meanRT_correct_left), meanRT_correct_right);
meanRT_incorrect=horzcat(fliplr(meanRT_incorrect_left), meanRT_incorrect_right);
SEMRT_correct=horzcat(fliplr(SEMRT_correct_left), SEMRT_correct_right);
SEMRT_incorrect=horzcat(fliplr(SEMRT_incorrect_left), SEMRT_incorrect_right);

meanslope_correct=horzcat(fliplr(meanslope_correct_left),meanslope_correct_right);
meanslope_incorrect=horzcat(fliplr(meanslope_incorrect_left), meanslope_incorrect_right);
SEMslope_correct=horzcat(fliplr(SEMslope_correct_left),SEMslope_correct_right);
SEMslope_incorrect=horzcat(fliplr(SEMslope_incorrect_left),SEMslope_incorrect_right);

meanacttime_correct=horzcat(fliplr(meanacttime_correct_left), meanacttime_correct_right);
meanacttime_incorrect=horzcat(fliplr(meanacttime_incorrect_left), meanacttime_incorrect_right);
SEMacttime_correct=horzcat(fliplr(SEMacttime_correct_left),SEMacttime_correct_right);
SEMacttime_incorrect=horzcat(fliplr(SEMacttime_incorrect_left),SEMacttime_incorrect_right);

xright=[c4 c3 c2 c1]; xleft=[1/c1 1/c2 1/c3 1/c4];
x=horzcat(xleft,xright);
xstr=sprintf('R/L %s ratio',taskstr);
figure(2)
subaxis(2,2,1,'SpacingVert',0.02,'MR',0.01); 
hold on
pc_control=plot(x,meanRT_correct(1,:),'-blueo','MarkerFaceColor','blue')
pi_control=plot(x,meanRT_incorrect(1,:),'-redo','MarkerFaceColor','red')
errorbar(x,meanRT_correct(1,:), SEMRT_correct(1,:), 'blueo')
errorbar(x,meanRT_incorrect(1,:), SEMRT_incorrect(1,:), 'redo')
xlim([0 102]); ylim([0 8000])
ylabel('Milliseconds after stim on'); 
ax = gca; ax.XScale = 'log'
set(gca,'XTickLabel', [])
set(gca,'box','off')
legend('Control correct', 'Control incorrect')
titlestr=sprintf('Average subjective reaction time by R/L %s ratio' , taskstr);
title(titlestr);

subaxis(2,2,3,'SpacingVert',0.02,'MR',0.01); 
hold on
pc_led=plot(x,meanRT_correct(2,:),'-cyano','MarkerFaceColor','cyan')
pi_led=plot(x,meanRT_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
errorbar(x,meanRT_correct(2,:), SEMRT_correct(2,:), 'cyano')
errorbar(x,meanRT_incorrect(2,:), SEMRT_incorrect(2,:), 'magentao')
xlabel(xstr); ylabel('Milliseconds after stim on'); 
xlim([0 102]);ylim([0 8000])
ax = gca; ax.XScale = 'log'
set(gca,'box','off')
legend('LED correct', 'LED incorrect')

subaxis(2,2,2,'SpacingVert',0.02,'MR',0.01); 
hold on
pcc_acttime=plot(x,meanacttime_correct(1,:),'-blueo','MarkerFaceColor','blue')
pci_acttime=plot(x,meanacttime_incorrect(1,:),'-redo','MarkerFaceColor','red')
errorbar(x,meanacttime_correct(1,:),SEMacttime_correct(1,:),'blueo','LineStyle','none');
errorbar(x,meanacttime_incorrect(1,:), SEMacttime_incorrect(1,:),'redo','LineStyle','none');
xlim([0 102]);ylim([0 8000])
ax = gca; ax.XScale = 'log'
set(gca,'XTickLabel',[],'YTickLabel', [])
set(gca,'box','off')
titlestr=sprintf('Average true decision time by R/L %s ratio' , taskstr);
title(titlestr);

subaxis(2,2,4,'SpacingVert',0.02,'MR',0.01); 
hold on
pledi_acttime=plot(x,meanacttime_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
pledc_acttime=plot(x,meanacttime_correct(2,:),'-cyano','MarkerFaceColor','cyan')
errorbar(x,meanacttime_correct(2,:),SEMacttime_correct(2,:),'cyano','LineStyle','none');
errorbar(x,meanacttime_incorrect(2,:), SEMacttime_incorrect(2,:),'magentao','LineStyle','none');
xlim([0 102]);ylim([0 8000])
xlabel(xstr);  
set(gca,'YTickLabel', [])
ax = gca; ax.XScale = 'log'
set(gca,'box','off')

figure(3)
subaxis(1,2,1,'SpacingHor',0.02,'MR',0.01); 
hold on
xlim([0 102]); ylim([0 0.2])
ax = gca
ax.XScale = 'log'
set(gca,'box','off')
xlabel(xstr); ylabel('Absolute value of quadrature trace slope'); 
titlestr=sprintf('Average quadrature response trace slope by R/L %s ratio.', taskstr);
title(titlestr)
pcc_rxnslope=plot(x,meanslope_correct(1,:),'-blueo','MarkerFaceColor','blue')
pci_rxnslope=plot(x,meanslope_incorrect(1,:),'-redo','MarkerFaceColor','red')
errorbar(x,meanslope_correct(1,:),SEMslope_correct(1,:),'blueo','LineStyle','none');
errorbar(x,meanslope_incorrect(1,:), SEMslope_incorrect(1,:),'redo','LineStyle','none');
legend('Control Correct', 'Control Incorrect')

subaxis(1,2,2,'SpacingHor',0.02,'MR',0.01); 
hold on
xlim([0 102]); ylim([0 0.2])
ax=gca
ax.XScale = 'log'
xlabel(xstr);
set(gca, 'YTickLabel',[])
set(gca,'box','off')
titlestr=sprintf('Average quadrature response trace slope by R/L %s ratio.' , taskstr);
title(titlestr)
pledc_rxnslope=plot(x,meanslope_correct(2,:),'-cyano','MarkerFaceColor','cyan')
pledi_rxnslope=plot(x,meanslope_incorrect(2,:),'-magentao','MarkerFaceColor','magenta')
errorbar(x,meanslope_correct(2,:),SEMslope_correct(2,:),'cyano','LineStyle','none');
errorbar(x,meanslope_incorrect(2,:), SEMslope_incorrect(2,:),'magentao','LineStyle','none');
legend('LED Correct', 'LED Incorrect')

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

               
               
               
               