%sort trials based on target contrast as well
%need to loop through trials of each contrast ratio for each target contrast first and get all 
%mean qvals
%find rxntime by side && contrast ratio with qvals_adj
%x axis vector will be 1/c1 1/c2 1/c3 1/c4 c4 c3 c2 c1
%run through 4 times to plot all found rxntimes per side and contrast ratio for each target contrast

%%
cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
savedate=input.saveTime;
ndatarun=1;
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
SIx = double(strcmp(input.trialOutcomeCell,'success'));
FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
block2=celleqel2mat_padded(input.tBlock2TrialNumber);
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg);
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity; rdt=-1*ldt;
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
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} (input.quadratureValues{trN}(end)+input.quadratureValues{trN+1})]);
        stimTime = double(input.stimTimestampMs{trN});
        qTimes_zero = qTimes-stimTime;
        qVals = qVals-qVals(1);
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000);
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind);
            qVals_sub = qVals(time_ind);
            qTimes_temp = qTimes(time_ind);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000;
            qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)';
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)';
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
ntgts=length(unqtargets);
%control=1 led=2, incorrect=1 correct=2, right=1 left=2, irat=1=c4
%irat=2=c3 irat=3=c2 irat=4=c1
sorted_trials = nan(nledcond, ncorrectcond, nstimlocation, ntgts, nratios, numtrials);
for icontrol=1:nledcond
  for ibehavior = 1:ncorrectcond
    for iloc = 1:nstimlocation
        for itgt=1:ntgts
            for irat = 1:nratios
               clear tempcon tempbehavesuccess tempbehavefailure tempbehave temploc temprat temp temptgt
               tempcon = find(block2 == icontrol-1);
               if ibehavior == 1
                   tempbehave = find(FIx == 1);
               elseif ibehavior == 2
                   tempbehave = find(SIx == 1);
               end
               temploc = find(tLeftTrial == iloc-1);
               temprat = find(contrastratio == uncontrasts(irat));
               temptgt = find(tContrast==unqtargets(itgt));
               temp =intersect(intersect(intersect(intersect(tempcon, tempbehave), temploc), temptgt),temprat);
               sorted_trials(icontrol, ibehavior, iloc, irat, 1:length(temp)) = temp;
           end; end; end; end; end;
%%
rxntime=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
for icontrol=1:nledcond
  for ibehavior = 1:ncorrectcond
    for iloc = 1:nstimlocation
        for itgt=1:ntgts
            for irat = 1:nratios
                clear temp temptrials temptime qvals_temp_mean
                 temp = squeeze(sorted_trials(icontrol,ibehavior,iloc,itgt,irat,:));
                 temptrials=temp(find(~isnan(temp)));
                 qvals_adj=qVals_final(:,temptrials)-qVals_final(8000,temptrials);
                 qvals_temp_mean=nanmean(qvals_adj,2);
                  for itrial=1:length(temptrials)
                      if iloc==1 && ibehavior==2 %right correct
                          temptime=min(find(qvals_temp_mean(8000:18000)<=rdt));
                      elseif iloc==1 && ibehavior==1 %right incorrect
                          temptime=min(find(qvals_temp_mean(8000:18000)>=ldt));
                      elseif iloc==2 && ibehavior==2 %left correct
                          temptime=min(find(qvals_temp_mean(8000:18000)>=ldt));
                      elseif iloc==2 && ibehavior==1 %left incorrect
                          temptime=min(find(qvals_temp_mean(8000:18000)<=rdt));
                      end    
                      if ~isempty(temptime)
                        rxntime(icontrol,ibehavior,iloc,itgt,irat,1)=temptime;  
                      end                               
            end; end;end;end;end; end

%% maybe?
x=[1/c1 1/c2 1/c3 1/c4 c4 c3 c2 c1];
unqtargets_label=round(unqtargets,2);
for itgt=1:ntgts
    for icontrol=1:nledcond
       for ibehavior=1:ncorrectcond
            for iloc=1:nstimlocation
                for irat=1:nratios
                    
                    if   ibehavior==2 
                        RTconcat_ctrl(1,ibehavior,itgt,1:4)=squeeze(rxntime(1,ibehavior,2,itgt,1:4));
                        RTconcat_ctrl(1,ibehavior,itgt,5:8)=squeeze(rxntime(1,ibehavior,1,itgt,1:4));
                        RTconcat_led(1,ibehavior,itgt,1:4)=squeeze(rxntime(2,ibehavior,2,itgt,1:4));
                        RTconcat_led(1,ibehavior,itgt,5:8)=squeeze(rxntime(2,ibehavior,1,itgt,1:4));
                        figure(1)
                        subaxis(2,2,itgt,'SpacingVert',0.02,'SpacingHor',0.02,'MR',0.01); 
                         hold on
                         titlestr=sprintf('Tcon= %g',unqtargets_label(itgt));
                         title(titlestr) 
                                                
                        r1=plot(x, squeeze(RTconcat_ctrl(1,ibehavior,itgt,:)),'-blueo', 'MarkerFaceColor','blue')
                        r2=plot(x,squeeze(RTconcat_led(1,ibehavior,itgt,:)),'-cyano','MarkerFaceColor','cyan')
                        xlim([0 102]); ylim([0 8000])
                        ax = gca; ax.XScale = 'log'; set(gca,'box','off')                       
                        if itgt==1 | itgt==2
                            set(gca, 'XTickLabel',[]);end                    
                        if itgt==2 || itgt==4
                            set(gca,'YTickLabel',[]);end                        
                        if itgt==1 |itgt==3                            
                            ylabel('Milliseconds');end                        
                        if itgt==3 | itgt==4
                            xlabel('R/L Contrast Ratio');end
                        if itgt==1
                            legend('Control correct', 'LED correct'); end
                    elseif ibehavior==1
                        figure(2)
                         RTconcat_ctrl_inc(1,ibehavior,itgt,1:4)=squeeze(rxntime(1,ibehavior,2,itgt,1:4));
                        RTconcat_ctrl_inc(1,ibehavior,itgt,5:8)=squeeze(rxntime(1,ibehavior,1,itgt,1:4));
                        RTconcat_led_inc(1,ibehavior,itgt,1:4)=squeeze(rxntime(2,ibehavior,2,itgt,1:4));
                        RTconcat_led_inc(1,ibehavior,itgt,5:8)=squeeze(rxntime(2,ibehavior,1,itgt,1:4));
                        subaxis(2,2,itgt,'SpacingVert',0.02,'SpacingHor', 0.02,'MR',0.01); 
                         hold on
                         titlestr=sprintf('Tcon= %g',unqtargets_label(itgt));
                         title(titlestr) 
                                                    
                        r1=plot(x, squeeze(RTconcat_ctrl_inc(1,ibehavior,itgt,:)),'-redo', 'MarkerFaceColor','red')
                        r2=plot(x,squeeze(RTconcat_led_inc(1,ibehavior,itgt,:)),'-magentao','MarkerFaceColor','magenta')
                        xlim([0 102]);ylim([0 8000])
                        ax = gca; ax.XScale = 'log'; set(gca,'box','off')
                        if itgt==1 | itgt==2
                            set(gca, 'XTickLabel',[]);end                     
                        if itgt==2 || itgt==4
                            set(gca,'YTickLabel',[]);end                        
                        if itgt==1 |itgt==3                            
                            ylabel('Milliseconds');end                        
                        if itgt==3 | itgt==4
                            xlabel('R/L Contrast Ratio');end                       
                        if itgt==1                            
                            legend('Control incorrect', 'LED incorrect'); end;                                        
                    end           
                end; end; end; end;  end


         
                        
                        
                     
 
                 

        
                   