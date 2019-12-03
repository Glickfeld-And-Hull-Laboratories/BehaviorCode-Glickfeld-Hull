%Run find_mult_files_for EGB_edited first to load all days
%Next run load_contrast_data_to_struct_forEGB_edited_602 to load all data
%to single structure 

%% For Single day of data (input) 
for i = 1:length(input)
Asuccesses = double(strcmp(input.trialOutcomeCell,'success'));
end
Asuccesses=Asuccesses(qVals_idx);

for i = 1:length(input)
Aincorrects = double(strcmp(input.trialOutcomeCell,'incorrect'));
end
Aincorrects=Aincorrects(qVals_idx);

for i = 1:length(input)
Aignores = double(strcmp(input.trialOutcomeCell,'ignore'));
end
Aignores=Aignores(qVals_idx);

Aori = double(cell2mat(input.tGratingDirectionStart));
Aori = Aori(qVals_idx);

Alefttrials = double(cell2mat(input.tGratingDirectionStart));
Alefttrials = Alefttrials(qVals_idx);

NoAdapt = double(cell2mat(input.aGratingContrast));
NoAdapt = NoAdapt(qVals_idx);
NoAdapt = ~NoAdapt;

Aadaptori = double(cell2mat(input.aGratingDirectionDeg)); 
Aadaptori = Aadaptori(qVals_idx);

A90Trials = Asuccesses(Aadaptori==90 & NoAdapt==0);
% A90Trials =A90Trials(A90_qVals_idx);
A90Oris = Aori(Aadaptori==90 & NoAdapt==0);
% A90Oris = A90Oris(A90_qVals_idx);
A90Ignores = Aignores(Aadaptori==90 & NoAdapt==0);
numIgnores90 = sum(A90Ignores);

A0Trials = Asuccesses(Aadaptori==0 & NoAdapt==0);
% A0Trials = A0Trials(A0_qVals_idx);
A0Oris = Aori(Aadaptori==0 & NoAdapt==0);
% A0Oris = A0Oris(A0_qVals_idx);
A0Ignores = (Aignores(Aadaptori==0 & NoAdapt==0));
numIgnores0 = sum(A0Ignores);

AnaTrials = Asuccesses(NoAdapt==1);
% AnaTrials = AnaTrials(Ana_qVals_idx);
AnaOris = Aori(NoAdapt==1);
% AnaOris = AnaOris(Ana_qVals_idx);
AnaIgnores = Aignores(NoAdapt==1);
numIgnoresNA = sum(AnaIgnores);

Auniqueori = unique(Aori);

A90success=nan(1,length(Auniqueori)); %A90success is array of # of total trials correct for each orientation??
A0success=nan(1,length(Auniqueori));
Anasuccess=nan(1,length(Auniqueori));


for Anori=1:length(Auniqueori)
    A90success(Anori) = nansum(A90Trials(A90Oris==Auniqueori(Anori)));
    NumIgnores90(Anori) = sum(Aignores(Aadaptori==90 & NoAdapt==0& Aori==Auniqueori(Anori)));

    A90alloritrials(Anori) = length(A90Trials(A90Oris==Auniqueori(Anori)))-NumIgnores90(Anori);
    if Auniqueori(Anori)>0
        A90pctL(Anori) = A90success(Anori)./A90alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A90pctL(Anori) = 1-(A90success(Anori)./A90alloritrials(Anori));
    end
end

%adding %Correct & 95%CIs for each ori in 0deg adapt condition
for Anori=1:length(Auniqueori)
    [A90pctCorrect(Anori),A90CI(Anori,:)]=binofit(A90success(Anori),A90alloritrials(Anori));
end
A90CI = A90CI';
A90CI = abs(A90pctCorrect-A90CI);

for Anori=1:length(Auniqueori)
    A0success(Anori) = nansum(A0Trials(A0Oris==Auniqueori(Anori)));
    NumIgnores0(Anori) = sum(Aignores(Aadaptori==0 & NoAdapt==0 & Aori==Auniqueori(Anori)));
    
    A0alloritrials(Anori) = length(A0Trials(A0Oris==Auniqueori(Anori)))-NumIgnores0(Anori);
    if Auniqueori(Anori)>0
        A0pctL(Anori) = A0success(Anori)./A0alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A0pctL(Anori) = 1-(A0success(Anori)./A0alloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [A0pctCorrect(Anori),A0CI(Anori,:)]=binofit(A0success(Anori),A0alloritrials(Anori)); 
end
A0CI = A0CI';
A0CI = abs(A0pctCorrect-A0CI);

for Anori=1:length(Auniqueori)
    Anasuccess(Anori) = nansum(AnaTrials(AnaOris==Auniqueori(Anori)));
    NumIgnoresNA(Anori) = sum(Aignores(NoAdapt==1 & Aori==Auniqueori(Anori)));
    
    Anaalloritrials(Anori) = length(AnaTrials(AnaOris==Auniqueori(Anori)))-NumIgnoresNA(Anori);
    if Auniqueori(Anori)>0
        AnapctL(Anori) = Anasuccess(Anori)./Anaalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        AnapctL(Anori) = 1-(Anasuccess(Anori)./Anaalloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [AnapctCorrect(Anori),AnaCI(Anori,:)]=binofit(Anasuccess(Anori),Anaalloritrials(Anori)); 
end
AnaCI = AnaCI';
AnaCI = abs(AnapctCorrect-AnaCI); 
 
figure; errorbar(Auniqueori,A90pctL,A90CI(1,:),A90CI(2,:),'r','LineWidth',1.5)
hold on
ylim([0,1])
errorbar(Auniqueori,A0pctL,A0CI(1,:),A0CI(2,:),'g','LineWidth',1.5)
errorbar(Auniqueori,AnapctL,AnaCI(1,:),AnaCI(2,:),'b','LineWidth',1.5)
title({'Psychometric Function per Adaptation Condition'; 'Paradigm: Flashing Adapt'; 'Subject: i1402'; 'n=6672 trials'})
xlabel('Target Orientation (degrees)')
ylabel('% Left Response')
legend('90 deg adapter, n=2241','0 deg adapter, n=2253','no adapter, n=2178')
hold off

% idx = logical([1 1 1 1 0 0 1 1 1 1]); 
% figure; errorbar(Auniqueori(idx),A90pctL(idx),A90CI(1,idx),A90CI(2,idx),'r','LineWidth',1.5)
% hold on
% ylim([0,1])
% errorbar(Auniqueori(idx),A0pctL(idx),A0CI(1,idx),A0CI(2,idx),'g','LineWidth',1.5)
% errorbar(Auniqueori(idx),AnapctL(idx),AnaCI(1,idx),AnaCI(2,idx),'b','LineWidth',1.5)
% title({'Psychometric Function per Adaptation Condition'; 'Paradigm: Flashing Adapt'; 'Subject: i1401'; 'n=4089 trials'})
% xlabel('Target Orientation (degrees)')
% ylabel('% Left Response')
% legend('90 deg adapter, n=1457','0 deg adapter, n=1496','no adapter, n=1137')
% hold off

%% For multiple days of data 
%First need to separate out adapt and non-adapt days 
doadapt = []; 
for i = 1:length(days_to_include)
    doadapt = [doadapt s(i).doAdapt];
end

%Next goal here to combine data from different days into single matrices. 
successes = [];
incorrects = [];
lefttrials = [];
ori = [];
Asuccesses = [];
Aincorrects = [];
Alefttrials = [];
Aadaptori = [];
Aori = [];
Acontrast = [];

for i = 1:length(days_to_include)
    if doadapt(i)==0
        successes = [successes s(i).SIx];
        incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
        lefttrials = [lefttrials s(i).tLeftTrial];
        ori = [ori s(i).tGratingDirectionStart];
    elseif doadapt(i)==1
        Asuccesses = [Asuccesses s(i).SIx];
        Aincorrects = [Aincorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
        Alefttrials = [Alefttrials s(i).tLeftTrial];
        Aadaptori = [Aadaptori s(i).aGratingDirectionDeg];
        Aori = [Aori s(i).tGratingDirectionStart];
        Acontrast = [Acontrast s(i).aGratingContrast];
    end
end



%Make matrix of ignores 
Aignores = [];
for i = 1:length(Asuccesses)
    if Asuccesses(i)==0 & Aincorrects(i)==0
        Aignores(i)=1;
    elseif Asuccesses(i)==1 | Aincorrects(i)==1
        Aignores(i)=0;
    end
end
AtotaltrialsNI = length(Asuccesses)-(sum(Aignores));

%throw out ignores (or rather, make nans) from "Asuccesses" and "Aincorrects" matrices - ANI =
%DoAdaptOn, No Ignores

% for i = 1:length(Asuccesses)
%     if Aignores(i)==1
%         ANIsuccesses(i)=nan;
%         ANIincorrects(i)=nan;
%     elseif Aignores(i)==0
%         ANIsuccesses(i)=Asuccesses(i);
%         ANIincorrects(i)=Aincorrects(i);
%     end
% end

%now have matrices of vars across all days, separated by adapt days and non-adapt days 

posoris = ori>0;
trialtypetest = sum(lefttrials-posoris);
Aposoris = Aori>0;
Atrialtypetest = sum(Alefttrials-Aposoris);
%+oris = Ltrials, -oris = Rtrials. Test that this aligns w tlefttrial var
%trialtypetest should =0


%%
%psych curve, percent Left 

%NoAdapt indicates trials in the no-adaptor condition, where DoAdapt was
%turned on, but the contrast was set to 0
NoAdapt=nan(1,length(Asuccesses));
for itrial=1:length(Asuccesses)
    if Acontrast(itrial)==1
        NoAdapt(itrial)=0;
    elseif Acontrast(itrial)==0
        NoAdapt(itrial)=1;
    end
end

%separate out for 0 and 90 degree adaptation conditions 
Auniqueori = unique(Aori);

%Don't include orientations with too-few trials (really small oris)
% for i=1:length(numOris)
%     if numOris(i)>-5 & numOris(i)<5
%         Auniqueori(i)=NaN;
%     elseif numOris(i)<-5 & numOris(i)>5
%         Auniqueori(i)=numOris(i);
%     end
% end
 
A90Trials = Asuccesses(Aadaptori==90 & NoAdapt==0);
A90Oris = Aori(Aadaptori==90 & NoAdapt==0);
numIgnores90 = sum(Aignores(Aadaptori==90 & NoAdapt==0));

A0Trials = Asuccesses(Aadaptori==0 & NoAdapt==0);
A0Oris = Aori(Aadaptori==0 & NoAdapt==0);
numIgnores0 = sum(Aignores(Aadaptori==0 & NoAdapt==0));

AnaTrials = Asuccesses(NoAdapt==1);
AnaOris = Aori(NoAdapt==1);
numIgnoresNA = sum(Aignores(NoAdapt==1));

A90success=nan(1,length(Auniqueori)); %A90success is array of # of total trials correct for each orientation??
A0success=nan(1,length(Auniqueori));
Anasuccess=nan(1,length(Auniqueori));

for Anori=1:length(Auniqueori)
    A90success(Anori) = nansum(A90Trials(A90Oris==Auniqueori(Anori)));
    NumIgnores90(Anori) = sum(Aignores(Aadaptori==90 & NoAdapt==0& Aori==Auniqueori(Anori)));

    A90alloritrials(Anori) = length(A90Trials(A90Oris==Auniqueori(Anori)))-NumIgnores90(Anori);
    if Auniqueori(Anori)>0
        A90pctL(Anori) = A90success(Anori)./A90alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A90pctL(Anori) = 1-(A90success(Anori)./A90alloritrials(Anori));
    end
end

%adding %Correct & 95%CIs for each ori in 0deg adapt condition
for Anori=1:length(Auniqueori)
    [A90pctCorrect(Anori),A90CI(Anori,:)]=binofit(A90success(Anori),A90alloritrials(Anori));
end
A90CI = A90CI';
A90CI = abs(A90pctCorrect-A90CI);

for Anori=1:length(Auniqueori)
    A0success(Anori) = nansum(A0Trials(A0Oris==Auniqueori(Anori)));
    NumIgnores0(Anori) = sum(Aignores(Aadaptori==0 & NoAdapt==0 & Aori==Auniqueori(Anori)));
    
    A0alloritrials(Anori) = length(A0Trials(A0Oris==Auniqueori(Anori)))-NumIgnores0(Anori);
    if Auniqueori(Anori)>0
        A0pctL(Anori) = A0success(Anori)./A0alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A0pctL(Anori) = 1-(A0success(Anori)./A0alloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [A0pctCorrect(Anori),A0CI(Anori,:)]=binofit(A0success(Anori),A0alloritrials(Anori)); 
end
A0CI = A0CI';
A0CI = abs(A0pctCorrect-A0CI);

for Anori=1:length(Auniqueori)
    Anasuccess(Anori) = nansum(AnaTrials(AnaOris==Auniqueori(Anori)));
    NumIgnoresNA(Anori) = sum(Aignores(NoAdapt==1 & Aori==Auniqueori(Anori)));
    
    Anaalloritrials(Anori) = length(AnaTrials(AnaOris==Auniqueori(Anori)))-NumIgnoresNA(Anori);
    if Auniqueori(Anori)>0
        AnapctL(Anori) = Anasuccess(Anori)./Anaalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        AnapctL(Anori) = 1-(Anasuccess(Anori)./Anaalloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    [AnapctCorrect(Anori),AnaCI(Anori,:)]=binofit(Anasuccess(Anori),Anaalloritrials(Anori)); 
end
AnaCI = AnaCI';
AnaCI = abs(AnapctCorrect-AnaCI); 
% 
% idx = logical([1 1 1 1 0 0 1 1 1 1]); 

figure; errorbar(Auniqueori,A90pctL,A90CI(1,:),A90CI(2,:),'r','LineWidth',1.5)
% figure; errorbar(Auniqueori(idx),A90pctL(idx),A90CI(1,idx),A90CI(2,idx),'r','LineWidth',1.5)
hold on
ylim([0,1])
% errorbar(Auniqueori(idx),A0pctL(idx),A0CI(1,idx),A0CI(2,idx),'g','LineWidth',1.5)
% errorbar(Auniqueori(idx),AnapctL(idx),AnaCI(1,idx),AnaCI(2,idx),'b','LineWidth',1.5)
errorbar(Auniqueori,A0pctL,A0CI(1,:),A0CI(2,:),'g','LineWidth',1.5)
errorbar(Auniqueori,AnapctL,AnaCI(1,:),AnaCI(2,:),'b','LineWidth',1.5)
title({'Psychometric Function per Adaptation Condition'; 'Paradigm: Flashing Adapt'; 'Subject: i1402'; 'n=5036 trials'})
xlabel('Target Orientation (degrees)')
ylabel('% Left Response')
legend('90 deg adapter, n=1764','0 deg adapter, n=1663','no adapter, n=1609')
hold off

%%
%plot ignores across conditions

numIgnores90 = sum(Aignores(Aadaptori==90 & NoAdapt==0));
numIgnores0 = sum(Aignores(Aadaptori==0 & NoAdapt==0));
numIgnoresNA = sum(Aignores(NoAdapt==1));

[ignores90Pct, ignores90CI] = binofit(numIgnores90,length(A90Trials));
[ignores0Pct, ignores0CI] = binofit(numIgnores0,length(A0Trials));
[ignoresNAPct, ignoresNACI] = binofit(numIgnoresNA,length(AnaTrials));

ignoresByTrialType = horzcat(ignoresNAPct,ignores90Pct,ignores0Pct);
ignores90CI = ignores90Pct-ignores90CI;
ignores0CI = ignores0Pct-ignores0CI;
ignoresNACI = ignoresNAPct-ignoresNACI;
ignoresCI = vertcat(ignoresNACI,ignores90CI,ignores0CI);
% x = ["90deg Adaptor","0deg Adaptor","No Adaptor"];

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});
figure
b=bar(conditions,ignoresByTrialType,'FaceColor',[1 1 1],'LineWidth',1.5)
hold on
b.EdgeColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
% set(gca,'XTick',1:3,'XTickLabel',x)
% ylim([0,0.03])
ylabel('ignore rate')
title({'Ignore Rate per Trial Type'; 'Flash Adapt'; 'Subject:i1402'})
errorbar(1:3,ignoresByTrialType,ignoresCI(:,1),ignoresCI(:,2),'k.','LineWidth',1.0)
hold off



%% Lapse rate/pct correct on easiest trials
AnaIgnores = Aignores(NoAdapt==1);
A90Ignores = Aignores(NoAdapt==0 & Aadaptori==90);
A0Ignores = Aignores(NoAdapt==0 & Aadaptori==0);

AnaEasyCorrects = sum(AnaTrials(AnaOris==45)) + sum(AnaTrials(AnaOris==-45));
AnaEasyTrials = (length(AnaTrials(AnaOris==45))-sum(AnaIgnores(AnaOris==45)))+ (length(AnaTrials(AnaOris==-45))-sum(AnaIgnores(AnaOris==-45)));
AnaEasyIncorrects = AnaEasyTrials-AnaEasyCorrects;
AnaPctEasyCorrect = AnaEasyCorrects/AnaEasyTrials;
AnaLapseRate = 1 - AnaPctEasyCorrect;

A90EasyCorrects = sum(A90Trials(A90Oris==45)) + sum(A90Trials(A90Oris==-45));
A90EasyTrials = (length(A90Trials(A90Oris==45))-sum(A90Ignores(A90Oris==45)))+ (length(A90Trials(A90Oris==-45))-sum(A90Ignores(A90Oris==-45)));
A90EasyIncorrects = A90EasyTrials-A90EasyCorrects;
A90PctEasyCorrect = A90EasyCorrects/A90EasyTrials;
A90LapseRate = 1 - A90PctEasyCorrect;

A0EasyCorrects = sum(A0Trials(A0Oris==45)) + sum(A0Trials(A0Oris==-45));
A0EasyTrials = (length(A0Trials(A0Oris==45))-sum(A0Ignores(A0Oris==45)))+(length(A0Trials(A0Oris==-45))-sum(A0Ignores(A0Oris==-45)));
A0EasyIncorrects = A0EasyTrials-A0EasyCorrects;
A0PctEasyCorrect = A0EasyCorrects/A0EasyTrials;
A0LapseRate = 1 - A0PctEasyCorrect;

[AnaLapseRate,AnaLapseCI] = binofit(AnaEasyIncorrects,AnaEasyTrials)
[A90LapseRate,A90LapseCI] = binofit(A90EasyIncorrects,A90EasyTrials)
[A0LapseRate,A0LapseCI] = binofit(A0EasyIncorrects,A0EasyTrials)

LapseRateByTrialType = horzcat(AnaLapseRate,A90LapseRate,A0LapseRate);

AnaLapseCI = AnaLapseRate-AnaLapseCI;
A90LapseCI = A90LapseRate-A90LapseCI;
A0LapseCI = A0LapseRate-A0LapseCI;

LapseCI = vertcat(AnaLapseCI, A90LapseCI, A0LapseCI);

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});
figure
b=bar(conditions,LapseRateByTrialType,'FaceColor',[1 1 1],'LineWidth',1.5)
hold on
b.EdgeColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
% set(gca,'XTick',1:3,'XTickLabel',x)
% ylim([0,0.03])
ylabel('Lapse Rate')
title({'Lapse Rate per Trial Type'; 'Flashing Adapt'; 'Subject:i1401'})
errorbar(1:3,LapseRateByTrialType,LapseCI(:,1),LapseCI(:,2),'k.','LineWidth',1.0)
hold off




%% Pct Correct

uniqueori = unique(ori);
success=nan(1,length(uniqueori));
pct=[];

for nori=1:length(uniqueori)
    success(nori)= sum(successes(ori==uniqueori(nori)));
    alloritrials(nori) = length(successes(ori==uniqueori(nori)));
    pct(nori) = success(nori)./alloritrials(nori);
end

%need to do same for adapt days, but separate out for 0 and 90 degree
%adaptation conditions 
Auniqueori = unique(Aori);
 
A90Trials = Asuccesses(Aadaptori==90 & NoAdapt==0);
A90Oris = Aori(Aadaptori==90 & NoAdapt==0);

A0Trials = Asuccesses(Aadaptori==0 & NoAdapt==0);
A0Oris = Aori(Aadaptori==0 & NoAdapt==0);

AnaTrials = Asuccesses(NoAdapt==1);
AnaOris = Aori(NoAdapt==1);

A90success=nan(1,length(Auniqueori));
A0success=nan(1,length(Auniqueori));
Anasuccess=nan(1,length(Auniqueori));

for Anori=1:length(Auniqueori)
    A90success(Anori) = nansum(A90Trials(A90Oris==Auniqueori(Anori)));
    A90alloritrials(Anori) = length(A90Trials(A90Oris==Auniqueori(Anori)));
    A90pct(Anori) = A90success(Anori)./A90alloritrials(Anori);
end

for Anori=1:length(Auniqueori)
    A0success(Anori) = nansum(A0Trials(A0Oris==Auniqueori(Anori)));
    A0alloritrials(Anori) = length(A0Trials(A0Oris==Auniqueori(Anori)));
    A0pct(Anori) = A0success(Anori)./A0alloritrials(Anori);
end

for Anori=1:length(Auniqueori)
    Anasuccess(Anori) = nansum(AnaTrials(AnaOris==Auniqueori(Anori)));
    Anaalloritrials(Anori) = length(AnaTrials(AnaOris==Auniqueori(Anori)));
    Anapct(Anori) = Anasuccess(Anori)./Anaalloritrials(Anori);
end

% figure; errorbar(Auniqueori(idx),A90pct(idx),A90CI(1,idx),A90CI(2,idx),'r','LineWidth',1.5)
% hold on
% ylim([0.5,1])
% errorbar(Auniqueori(idx),A0pct(idx),A0CI(1,idx),A90CI(2,idx),'g','LineWidth',1.5)
% errorbar(Auniqueori(idx),Anapct(idx),AnaCI(1,idx),AnaCI(2,idx),'b','LineWidth',1.5)
figure; errorbar(Auniqueori,A90pct,A90CI(1,:),A90CI(2,:),'r','LineWidth',1.5)
hold on
ylim([0.5,1])
errorbar(Auniqueori,A0pct,A0CI(1,:),A90CI(2,:),'g','LineWidth',1.5)
errorbar(Auniqueori,Anapct,AnaCI(1,:),AnaCI(2,:),'b','LineWidth',1.5)
title({'Pct Correct per Adaptation Condition'; 'Paradigm: Flashing Adapt'; 'Subjects: i1402'; 'n=3856 trials'})
xlabel('Target Orientation (degrees)')
ylabel('% Correct Response')
legend('90 deg adapter, n=1305','0 deg adapter, n=1301','no adapter, n=1250')
hold off
hold off


