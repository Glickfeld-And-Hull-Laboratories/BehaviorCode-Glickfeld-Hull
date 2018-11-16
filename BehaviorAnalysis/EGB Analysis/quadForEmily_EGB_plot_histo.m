%Drag single day data file over

Iix = find(strcmp(input.trialOutcomeCell, 'ignore')); %command strcmp says match "trialOutcomeCell" when it says "ignore"
%Find all ignore trials

%Nix = find(celleqel2mat_padded(input.didNoGo)); 

Tix = setdiff(1:length(input.trialOutcomeCell), Iix); %"setdiff" means give 2 inputs, what is in 1 but not the other? 
%Creating an array of not-ignores

LeftTrials = double(cell2mat(input.tLeftTrial));
TLeftTrials = LeftTrials(find(Tix));
LeftTrials = TLeftTrials';

%%
%ADDING ABILITY TO DIVIDE BY SUCCESS AND INCORRECT
TrialOutcome = input.trialOutcomeCell;
TTrialOutcome = TrialOutcome(find(Tix));
nTrials = length(Tix);
TrialOutcomeMat = nan(nTrials,1);
for itrial = 1:nTrials
    clear temp 
    temp = TTrialOutcome{itrial};
    if strcmp(temp, 'success')  %'' means string; ==1 because if loop, has to be condition satisfying
        TrialOutcomeMat(itrial) = 1;
    elseif strcmp(temp, 'incorrect') 
        TrialOutcomeMat(itrial) = 0;
    end 
end

%maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
%Max decision time; will be single value 

qVals_final = nan(18001, uint16(length(TTrialOutcome))); %nan creates matrix without numbers, good for pre-allocating
%6 sec ITI, 2 sec stationary, 10 sec to decide- made cell for every ms +1, 18001
%uint16 is 16-bit unsigned integer 
%Essentially preparing or pre-allocatng matrix of nans for final quad value
%for each trial 

qTimes_act = nan(18001, uint16(length(TTrialOutcome)));
%act = actual, Qtime when QVals recorded
%have to "INTERPOLATE" qTime onto Mathworks time , mapping discrete times on
%continuous 
%qTimes_act prepares matrix for storing actual qTimes later

qTimes_thresh = nan(1, uint16(length(TTrialOutcome)));
%unknown

%cVals_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
%unknown

%%
%trN is "trial number" variable created for for:loop

for trN = 1:input.trialSinceReset-1 %saves need to pre-allocate because just builds as you go
    
    if find(Tix == trN) %if trial number you're on in that loop is not an ignore, then
        
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        %Build 1 row array in a matrix (builds row with each loop)
        %{} is for referencing numbers in a CELL vs matrix 
        %second column within row will be next trial
        %[] for putting multiple values into an array
        %/1000 for ms not us
        
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        %same thing for quad values
        
        %DONT NEED. cTimes = double(input.counterTimesUs{trN}./1000);
        %DONT NEED. cVals = double(input.counterValues{trN});
        
        stimTime = double(input.stimTimestampMs{trN});
        %read in mWorks time that stim comes on
        
        %DONT NEED. stimVal = double(input.qStimOn{trN});
        
        qTimes_zero = qTimes-stimTime;
        %For whole rows, 0 to stimTime on
        
        qVals = qVals-qVals(1);
        clear x
        x=find(qVals>1000);
        for i=1:length(x)
            qVals(x(i)) = (qVals(x(i)+1)+qVals(x(i)-1))/2;
        end
        clear x
        
        x=find(qVals<-1000);
        for i=1:length(x)
            qVals(x(i)) = (qVals(x(i)+1)+qVals(x(i)-1))/2;
        end
        
        %Zeroing to where qVal started on that trial 
        
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000);
        %Limiting qTimes to mWorks 180000 time limit - only interested in
        %this time (not sure what other time there is)
        
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
               if isempty(qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN))
                    qTimes_thresh(:,trN) = NaN;
                else
                    qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN);
                end
                %cVals_thresh(:,trN) = cVals(:,find(cTimes>qTimes_thresh(:,trN),1,'first'));
            end
        else
            return
        end
    end
end

qVals_final = qVals_final';

%%Main outputs are qVals_final - matrix of 180000 ms, each trial 1 column,
%%gives qValue for each trial at each ms

%Other main output is qTimes_final -8000 --> 10000; need to be able to plot
%or map values to certain time 

%Each trial is a single column in qVals, QTimes is just one long row

%QTimes, 0 is when stim comes on; shouldn't be hard to find Q value sright
%before, store and plot - think about which rows/cells to target for
%plotting

%%                               
CorrLTrials = find(LeftTrials==1 & TrialOutcomeMat==1);%output is trial #s that satisfy these conditions 
IncorrLTrials = find(LeftTrials==1 & TrialOutcomeMat==0);
CorrRTrials = find(LeftTrials==0 & TrialOutcomeMat==1);
IncorrRTrials = find(LeftTrials==0 & TrialOutcomeMat==0);

AllCorrTrials = find(TrialOutcomeMat==1);
AllIncorrTrials = find(TrialOutcomeMat==0);
AllTrials = find(TrialOutcomeMat==1|0);

qVals_CorrLTrials = qVals_final(CorrLTrials,:);
qVals_IncorrLTrials = qVals_final(IncorrLTrials,:);
qVals_CorrRTrials = qVals_final(CorrRTrials,:);
qVals_IncorrRTrials = qVals_final(IncorrRTrials,:);

qVals_AllCorr = qVals_final(AllCorrTrials,:);
qVals_AllIncorr = qVals_final(AllIncorrTrials,:);
qVals_AllTrials = qVals_final(AllTrials,:);

qVals_CorrL_4sec = qVals_CorrLTrials(:,find(qTimes_final==-4000):find(qTimes_final==0));
qVals_IncorrL_4sec = qVals_IncorrLTrials(:,find(qTimes_final==-4000):find(qTimes_final==0));
qVals_CorrR_4sec = qVals_CorrRTrials(:,find(qTimes_final==-4000):find(qTimes_final==0));
qVals_IncorrR_4sec = qVals_IncorrRTrials(:,find(qTimes_final==-4000):find(qTimes_final==0));

qVals_AllCorr_4sec = qVals_AllCorr(:,find(qTimes_final==-4000):find(qTimes_final==0));
qVals_AllIncorr_4sec = qVals_AllIncorr(:,find(qTimes_final==-4000):find(qTimes_final==0));
qVals_AllTrials_4sec = qVals_AllTrials(:,find(qTimes_final==-4000):find(qTimes_final==0));

%figure; plot(qVals_AllCorr_2sec(1,:))
%figure; plot(diff(qVals_AllCorr_2sec(1,:)))

%QdiffTest = diff(qVals_AllCorr_2sec,1,2);
%QdiffSum = sum(abs(QdiffTest),2);
%%
qVals_CorrL_sum = sum(abs(diff(qVals_CorrL_4sec,1,2)),2); 
qVals_IncorrL_sum = sum(abs(diff(qVals_IncorrL_4sec,1,2)),2);
qVals_CorrR_sum = sum(abs(diff(qVals_CorrR_4sec,1,2)),2);
qVals_IncorrR_sum = sum(abs(diff(qVals_IncorrR_4sec,1,2)),2);

qVals_AllCorr_sum = sum(abs(diff(qVals_AllCorr_4sec,1,2)),2);
qVals_AllIncorr_sum = sum(abs(diff(qVals_AllIncorr_4sec,1,2)),2);
qVals_AllTrials_sum = sum(abs(diff(qVals_AllTrials_4sec,1,2)),2);

%%
figure(1)
cdfplot(qVals_CorrL_sum)

figure(2)
cdfplot(qVals_IncorrL_sum)

figure(3)
cdfplot(qVals_CorrR_sum)

figure(4)
cdfplot(qVals_IncorrR_sum)

figure(5)
cdfplot(qVals_AllCorr_sum)

figure(6)
cdfplot(qVals_AllIncorr_sum)

figure(7)
cdfplot(qVals_AllTrials_sum)

%%
%plot histograms for correct/incorrect L and R trials, percentage of
%trials that fall into bins of absolute value sums of total movement during
%2 seconds prior to stim on time 

figure(1)
histogram(qVals_CorrL_sum, 'BinWidth',10, 'Normalization', 'probability')
xlim([0 200])
ylim([0 0.5])

figure(2)
histogram(qVals_IncorrL_sum, 'BinWidth', 10, 'Normalization','probability')
xlim([0 200])
ylim([0 0.5])

figure(3)
histogram(qVals_CorrR_sum, 'BinWidth', 10, 'Normalization', 'probability')
xlim([0 200])
ylim([0 0.5])

figure(4)
histogram(qVals_IncorrR_sum, 'BinWidth', 10, 'Normalization','probability')
xlim([0 200])
ylim([0 0.5])

figure(5)
histogram(qVals_AllCorr_sum, 'BinWidth', 10, 'Normalization', 'probability')
xlim([0 200])
ylim([0 0.5])

figure(6)
histogram(qVals_AllIncorr_sum, 'BinWidth', 10, 'Normalization', 'probability')
xlim([0 200])
ylim([0 0.5])

figure(7)
histogram(qVals_AllTrials_sum, 'BinWidth', 10, 'Normalization', 'probability')
xlim([0 200])
ylim([0 0.5])

%%
%TrueLeftTrials = find(LeftTrials==1);
%TrueRightTrials = find(LeftTrials==0);
%%
