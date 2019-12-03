%load input of vars & data, single day or collapsed 

for i = 1:length(input.trialOutcomeCell)
successes = double(strcmp(input.trialOutcomeCell,'success'));
end

clear i;

for i = 1:length(input.trialOutcomeCell)
ignores = double(strcmp(input.trialOutcomeCell,'ignore'));
end

clear i;

%trialOutcomeCell will be matrix of all trials with 1's as successes, 0's
%as incorrects, and 2's as ignores 
trialOutcomeCell = successes;
for i = 1:length(ignores)
    if ignores(i)==1
        trialOutcomeCell(i)=nan;
    elseif ignores(i)==0
        trialOutcomeCell(i)=trialOutcomeCell(i);
    end
end


PrevTrialOutcome = [nan trialOutcomeCell(1:4997)];

[tbl,chi2,p]=crosstab(trialOutcomeCell,PrevTrialOutcome)
tbl(1,3) = tbl(1,1) + tbl(1,2);
tbl(2,3) = tbl(2,1) + tbl(2,2);
tbl(3,1) = tbl(1,1) + tbl(2,1);
tbl(3,2) = tbl(1,2) + tbl(2,2);
tbl(3,3) = tbl(1,3) + tbl(2,3);

tblExpected = [(tbl(1,3)*tbl(3,1))/tbl(3,3) (tbl(1,3)*tbl(3,2))/tbl(3,3);...
    (tbl(2,3)*tbl(3,1))/tbl(3,3) (tbl(2,3)*tbl(3,2))/tbl(3,3)]

%Chi2 test
[tbl,chi2,p]=crosstab(trialOutcomeCell,PrevTrialOutcome)

%Pearson's correlation
X = trialOutcomeCell';
Y = PrevTrialOutcome';

Z = horzcat(Y,X);
R = corrcoef(Z,'Rows','complete');

%% Compare psych curves with diff outcomes on prev trials 

PrevTrialOutcome(isnan(PrevTrialOutcome))=2; %Prev trial outcome categories are 1=success, 0=incorrect, 2=ignore

Prev_Success_idx = find(PrevTrialOutcome==1);
Prev_Incorrect_idx = find(PrevTrialOutcome==0);
Prev_Ignore_idx = find(PrevTrialOutcome==2);

for i = 1:length(input)
Asuccesses = double(strcmp(input.trialOutcomeCell,'success'));
end
Asuccesses_PS=Asuccesses(Prev_Success_idx);
Asuccesses_PIN=Asuccesses(Prev_Incorrect_idx);
Asuccesses_PIG=Asuccesses(Prev_Ignore_idx);

for i = 1:length(input)
Aincorrects = double(strcmp(input.trialOutcomeCell,'incorrect'));
end
Aincorrects_PS=Aincorrects(Prev_Success_idx);
Aincorrects_PIN=Aincorrects(Prev_Incorrect_idx);
Aincorrects_PIG=Aincorrects(Prev_Ignore_idx);


for i = 1:length(input)
Aignores = double(strcmp(input.trialOutcomeCell,'ignore'));
end
Aignores_PS=Aignores(Prev_Success_idx);
Aignores_PIN=Aignores(Prev_Incorrect_idx);
Aignores_PIG=Aignores(Prev_Ignore_idx);


Aori = double(cell2mat(input.tGratingDirectionStart));
Aori_PS=Aori(Prev_Success_idx);
Aori_PIN=Aori(Prev_Incorrect_idx);
Aori_PIG=Aori(Prev_Ignore_idx);

Alefttrials = double(cell2mat(input.tGratingDirectionStart));
Alefttrials_PS=Alefttrials(Prev_Success_idx);
Alefttrials_PIN=Alefttrials(Prev_Incorrect_idx);
Alefttrials_PIG=Alefttrials(Prev_Ignore_idx);

NoAdapt = double(cell2mat(input.aGratingContrast));
NoAdapt_PS = NoAdapt(Prev_Success_idx);
NoAdapt_PIN = NoAdapt(Prev_Incorrect_idx);
NoAdapt_PIG = NoAdapt(Prev_Ignore_idx);
NoAdapt_PS = ~NoAdapt_PS;
NoAdapt_PIN = ~NoAdapt_PIN;
NoAdapt_PIG = ~NoAdapt_PIG;

Aadaptori = double(cell2mat(input.aGratingDirectionDeg)); 
Aadaptori_PS = Aadaptori(Prev_Success_idx);
Aadaptori_PIN = Aadaptori(Prev_Incorrect_idx);
Aadaptori_PIG = Aadaptori(Prev_Ignore_idx);

numIgnores_PS = sum(Aignores_PS);
numIgnores_PIN = sum(Aignores_PIN);
numIgnores_PIG = sum(Aignores_PIG);

Auniqueori = unique(Aori);

PS_success=nan(1,length(Auniqueori)); %PS_success is array of # of total trials correct for each orientation, in the "prev success" condition
PIN_success=nan(1,length(Auniqueori));
PIG_success=nan(1,length(Auniqueori));


for Anori=1:length(Auniqueori)
    PS_success(Anori) = nansum(Asuccesses_PS(Aori_PS==Auniqueori(Anori)));
    NumIgnoresPS(Anori) = sum(Aignores_PS(Aori_PS==Auniqueori(Anori)));
    APSalloritrials(Anori) = length(Asuccesses_PS(Aori_PS==Auniqueori(Anori)))-NumIgnoresPS(Anori);
    if Auniqueori(Anori)>0
        APSpctL(Anori) = PS_success(Anori)./APSalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        APSpctL(Anori) = 1-(PS_success(Anori)./APSalloritrials(Anori));
    end
end

%adding %Correct & 95%CIs for each ori in 0deg adapt condition
for Anori=1:length(Auniqueori)
    [APSpctCorrect(Anori),APSCI(Anori,:)]=binofit(PS_success(Anori),APSalloritrials(Anori));
end
APSCI = APSCI';
APSCI = abs(APSpctCorrect-APSCI);

for Anori=1:length(Auniqueori)
    PIN_success(Anori) = nansum(Asuccesses_PIN(Aori_PIN==Auniqueori(Anori)));
    NumIgnoresPIN(Anori) = sum(Aignores_PIN(Aori_PIN==Auniqueori(Anori)));
    APINalloritrials(Anori) = length(Asuccesses_PIN(Aori_PIN==Auniqueori(Anori)))-NumIgnoresPIN(Anori);
    if Auniqueori(Anori)>0
        APINpctL(Anori) = PIN_success(Anori)./APINalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        APINpctL(Anori) = 1-(PIN_success(Anori)./APINalloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [APINpctCorrect(Anori),APINCI(Anori,:)]=binofit(PIN_success(Anori),APINalloritrials(Anori)); 
end
APINCI = APINCI';
APINCI = abs(APINpctCorrect-APINCI);

for Anori=1:length(Auniqueori)
    PIG_success(Anori) = nansum(Asuccesses_PIG(Aori_PIG==Auniqueori(Anori)));
    NumIgnoresPIG(Anori) = sum(Aignores_PIG(Aori_PIG==Auniqueori(Anori)));
    APIGalloritrials(Anori) = length(Asuccesses_PIG(Aori_PIG==Auniqueori(Anori)))-NumIgnoresPIG(Anori);
    if Auniqueori(Anori)>0
        APIGpctL(Anori) = PIG_success(Anori)./APIGalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        APIGpctL(Anori) = 1-(PIG_success(Anori)./APIGalloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [APIGpctCorrect(Anori),APIGCI(Anori,:)]=binofit(PIG_success(Anori),APIGalloritrials(Anori)); 
end
APIGCI = APIGCI';
APIGCI = abs(APIGpctCorrect-APIGCI); 
 

figure; errorbar(Auniqueori,APSpctL,APSCI(1,:),APSCI(2,:),'r','LineWidth',1.5)
hold on
ylim([0,1])
errorbar(Auniqueori,APINpctL,APINCI(1,:),APINCI(2,:),'g','LineWidth',1.5)
errorbar(Auniqueori,APIGpctL,APIGCI(1,:),APIGCI(2,:),'b','LineWidth',1.5)
title({'Psychometric Functions per Prev Trial Outcome';'Subject: i1402'; 'n=4893 trials'})
xlabel('Target Orientation (degrees)')
ylabel('% Left Response')
legend('Previous Success, n=4195','Previous Incorrect, n=632','Previous Ignore, n=66')
hold off





