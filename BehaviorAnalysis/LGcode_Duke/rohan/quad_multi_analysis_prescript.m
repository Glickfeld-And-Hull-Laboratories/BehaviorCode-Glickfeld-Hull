cd('C:\Users\rohan\Documents\MATLAB\Repositories\BehaviorCode-Glickfeld-Hull-Master\BehaviorAnalysis')
ndatarun=1; %NEED TO CHANGE THIS FOR EACH RUN TO STORE VALUES ACROSS DAYS
Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
savedate=input.saveTime;
qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
SIx = double(strcmp(input.trialOutcomeCell,'success'));
FIx = double(strcmp(input.trialOutcomeCell,'incorrect'));
block2=celleqel2mat_padded(input.tBlock2TrialNumber);
eccentricity=celleqel2mat_padded(input.tGratingEccentricityDeg);
ldt=abs(eccentricity(2))/input.feedbackMotionSensitivity; rdt=-1*ldt;
tContrast=celleqel2mat_padded(input.tGratingContrast); % comment out for size discrim
dContrast=celleqel2mat_padded(input.dGratingContrast); % comment out for size discrim
%tContrast=celleqel2mat_padded(input.tGratingDiameterDeg); 
%dContrast=celleqel2mat_padded(input.dGratingDiameterDeg);
unqtargets=unique(tContrast);
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
rxntime=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxntime_avg=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnqval=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
rxnslope=nan(nledcond,ncorrectcond,nstimlocation,nratios,1);
for icontrol=1:nledcond
   for ibehavior = 1:ncorrectcond
       for iloc = 1:nstimlocation
           for irat = 1:nratios
               clear temp qvals_adj qvals_temp_mean time_adj temptime
               temp = squeeze(sorted_trials(icontrol,ibehavior,iloc,irat,:));               
               temp=temp(find(~isnan(temp)));
               qvals_adj=qVals_final(:,temp)-qVals_final(8000,temp);               
               qvals_temp_mean=nanmean(qvals_adj,2);
                  for itrial=1:length(temp)
                      if iloc==1 && ibehavior==2 %right correct
                          temptime(itrial)=min(find(qvals_adj(8000:18000)<=rdt));
                      elseif iloc==1 && ibehavior==1 %right incorrect
                          temptime(itrial)=min(find(qvals_adj(8000:18000)>=ldt));
                      elseif iloc==2 && ibehavior==2 %left correct
                          temptime(itrial)=min(find(qvals_adj(8000:18000)>=ldt));
                      elseif iloc==2 && ibehavior==1 %left incorrect
                          temptime(itrial)=min(find(qvals_adj(8000:18000)<=rdt));
                      end                      
                   rxntime(icontrol,ibehavior,iloc,irat,1:length(temptime))=temptime;
                   rxntime_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(temptime);
                   rxnslopewindowstart=rxntime(icontrol,ibehavior,iloc,irat,:)-50;
                   rxnslopewindowend=rxntime(icontrol,ibehavior,iloc,irat,:)+400;
                   rxnqval(icontrol,ibehavior,iloc,irat,1:length(temptime))=qvals_temp_mean(temptime+8000);    
                   rxnqval_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnqval(icontrol,ibehavior,iloc,irat,:));
                   rxnslope(icontrol,ibehavior,iloc,irat,1:length(temptime))=abs((qvals_adj(rxnslopewindowend)-qvals_adj(rxnslopewindowstart)))/450;
                   rxnslope_avg(icontrol,ibehavior,iloc,irat,1)=nanmean(rxnslope(icontrol,ibehavior,iloc,irat,:));
       
                  end
               hold on
           end; end; end; end           
%%
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
           if iloc==1               
                rxntime_correct_right(icontrol,irat,ndatarun)=temptime_correct;
                rxntime_incorrect_right(icontrol,irat,ndatarun)=temptime_incorrect;  
                rxnslope_correct_right(icontrol,irat,ndatarun)=tempslope_correct;
                rxnslope_incorrect_right(icontrol,irat,ndatarun)=tempslope_incorrect;
                acttime_correct_right(icontrol,irat,ndatarun)=acttime_correct;
                acttime_incorrect_right(icontrol,irat,ndatarun)=acttime_incorrect;
           elseif iloc==2
                rxntime_correct_left(icontrol,irat,ndatarun)=temptime_correct;
                rxntime_incorrect_left(icontrol,irat,ndatarun)=temptime_incorrect;  
                rxnslope_correct_left(icontrol,irat,ndatarun)=tempslope_correct;
                rxnslope_incorrect_left(icontrol,irat,ndatarun)=tempslope_incorrect;
                acttime_correct_left(icontrol,irat,ndatarun)=acttime_correct;
                acttime_incorrect_left(icontrol,irat,ndatarun)=acttime_incorrect;
           end
        end; end; end;
%%
meanslope_correct_right=nanmean(rxnslope_correct_right,3);
meanslope_correct_left=nanmean(rxnslope_correct_left,3);
meanslope_incorrect_right=nanmean(rxnslope_incorrect_right,3);
meanslope_incorrect_left=nanmean(rxnslope_incorrect_left,3);
meanRT_correct_right=nanmean(rxntime_correct_right,3);
meanRT_correct_left=nanmean(rxntime_correct_left,3);
meanRT_incorrect_right=nanmean(rxntime_incorrect_right,3);
meanRT_incorrect_left=nanmean(rxntime_incorrect_left,3);
meanacttime_correct_right=nanmean(acttime_correct_right,3);
meanacttime_correct_left=nanmean(acttime_correct_left,3);
meanacttime_incorrect_right=nanmean(acttime_incorrect_right,3);
meanacttime_incorrect_left=nanmean(acttime_incorrect_left,3);

SEMslope_correct_right=nanstd(rxnslope_correct_right,[],3)/(sqrt(length(rxnslope_correct_right)));
SEMslope_correct_left=nanstd(rxnslope_correct_left,[],3)/(sqrt(length(rxnslope_correct_left)));
SEMslope_incorrect_right=nanstd(rxnslope_incorrect_right,[],3)/(sqrt(length(rxnslope_incorrect_right)));
SEMslope_incorrect_left=nanstd(rxnslope_incorrect_left,[],3)/(sqrt(length(rxnslope_incorrect_left)));
SEMRT_correct_right=nanstd(rxntime_correct_right,[],3)/(sqrt(length(rxntime_correct_right)));
SEMRT_correct_left=nanstd(rxntime_correct_left,[],3)/(sqrt(length(rxntime_correct_left)));
SEMRT_incorrect_right=nanstd(rxntime_incorrect_right,[],3)/(sqrt(length(rxntime_incorrect_right)));
SEMRT_incorrect_left=nanstd(rxntime_incorrect_left,[],3)/(sqrt(length(rxntime_incorrect_left)));
SEMacttime_correct_right=nanstd(acttime_correct_right,[],3)/(sqrt(length(acttime_correct_right)));
SEMacttime_correct_left=nanstd(acttime_correct_left,[],3)/(sqrt(length(acttime_correct_left)));
SEMacttime_incorrect_right=nanstd(acttime_incorrect_right,[],3)/(sqrt(length(acttime_incorrect_right)));
SEMacttime_incorrect_left=nanstd(acttime_incorrect_left,[],3)/(sqrt(length(acttime_incorrect_left)));
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
