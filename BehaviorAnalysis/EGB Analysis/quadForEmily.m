%Drag single day data file over

Iix = find(strcmp(input.trialOutcomeCell, 'ignore')); %command strcmp says match "trialOutcomeCell" when it says "ignore"
%Find all ignore trials

%Nix = find(celleqel2mat_padded(input.didNoGo)); 

Tix = setdiff(1:length(input.trialOutcomeCell), Iix); %"setdiff" means give 2 inputs, what is in 1 but not the other? 
%Creating an array of not-ignores

maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
%Max decision time; will be single value 

qVals_final = nan(18001, uint16(length(input.trialOutcomeCell))); %nan creates matrix without numbers, good for pre-allocating
%6 sec ITI, 2 sec stationary, 10 sec to decide- made cell for every ms +1, 18001
%uint16 is 16-bit unsigned integer 
%Essentially preparing or pre-allocatng matrix of nans for final quad value
%for each trial 

qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
%act = actual, Qtime when QVals recorded
%have to "INTERPOLATE" qTime onto Mathworks time , mapping discrete times on
%continuous 
%qTimes_act prepares matrix for storing actual qTimes later

qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
%unknown

%cVals_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
%unknown

%%

%trN is "trial number" variable created for for:loop

for trN = 1:length(input.trialOutcomeCell)-1 %saves need to pre-allocate because just builds as you go
    
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
                qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN);
                %cVals_thresh(:,trN) = cVals(:,find(cTimes>qTimes_thresh(:,trN),1,'first'));
            end
        else
            return
        end
    end
end

%%Main outputs are qVals_final - matrix of 180000 ms, each trial 1 column,
%%gives qValue for each trial at each ms

%Other main output is qTimes_final -8000 --> 10000; need to be able to plot
%or map values to certain time 

%Each trial is a single column in qVals, QTimes is just one long row

%QTimes, 0 is when stim comes on; shouldn't be hard to find Q value sright
%before, store and plot - think about which rows/cells to target for
%plotting

%STOP HERE AND USE ROHAN'S quadTestScript_singleday_RB to plot 
%%
minR = input.tooFastTimeMs;
figure;
for it = 1:30
    subplot(5,6,it)
    plot(qTimes_final, qVals_final(:,it))
    xlim([-100 input.tDecisionTimeMs{it}])
    vline(minR)
    if input.tLeftTrial{it}
        title(['Left ' input.trialOutcomeCell{it}])
    else
        title(['Right ' input.trialOutcomeCell{it}])
    end
end


SIx = setdiff(find(strcmp(input.trialOutcomeCell, 'success')),Nix);
MIx = setdiff(find(strcmp(input.trialOutcomeCell, 'incorrect')),Nix);
left = celleqel2mat_padded(input.tLeftTrial);
maxR = input.reactionTimeMs;
lowR = intersect(find(celleqel2mat_padded(input.tDecisionTimeMs)<maxR), find(celleqel2mat_padded(input.tDecisionTimeMs)>minR));

qVals_offset = bsxfun(@minus, qVals_final, qVals_final(8001,:));
figure;
subplot(2,2,1)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(SIx,find(left)))))
xlim([-500 2000])
vline(minR)
title('Correct Left Trials')

subplot(2,2,2)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(SIx,find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Correct Right Trials')

subplot(2,2,3)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(MIx,find(left)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Left Trials')

subplot(2,2,4)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(MIx,find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Right Trials')

suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])

% print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome.pdf']), '-dpdf')

figure
% con_str = strvcat('g', 'c', 'b', 'k');
% cons(find(cons<0.0001)) = [];
% ncon = length(cons);
subplot(2,2,1)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(SIx, find(left)))),2), nanstd(qVals_offset(:,(intersect(SIx, find(left)))),[],2)./sqrt(length(intersect(SIx, find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(SIx, find(left==0)))),2), nanstd(qVals_offset(:,(intersect(SIx, find(left==0)))),[],2)./sqrt(length(intersect(SIx, find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all correct trials')

subplot(2,2,2)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(MIx, find(left)))),2), nanstd(qVals_offset(:,(intersect(MIx, find(left)))),[],2)./sqrt(length(intersect(MIx, find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(MIx, find(left==0)))),2), nanstd(qVals_offset(:,(intersect(MIx, find(left==0)))),[],2)./sqrt(length(intersect(MIx, find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all incorrect trials')

% subplot(2,2,3)
% for icon = 1:ncon
%     ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(SIx, find(left))));
%     if length(ind)>0
%         shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
%         hold on;
%     end
% end
% xlim([-500 2000])
% vline(minR)
% title('Avg left correct trials by contrast')
% 
% subplot(2,2,4)
% for icon = 1:ncon
%     ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(SIx, find(left==0))));
%     shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
%     hold on
% end
% xlim([-500 2000])
% vline(minR)
% title('Avg right correct trials by contrast')
suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])
% print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome_avg.pdf']), '-dpdf')
