%% Ignore and lapse rates across animals

%Begin by loading vars from s structure
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

Auniqueori = unique(Aori);

 
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

%% Ignore rates
numIgnores90 = sum(Aignores(Aadaptori==90 & NoAdapt==0));
numIgnores0 = sum(Aignores(Aadaptori==0 & NoAdapt==0));
numIgnoresNA = sum(Aignores(NoAdapt==1));

[ignores90Pct, ignores90CI] = binofit(numIgnores90,length(A90Trials));
[ignores0Pct, ignores0CI] = binofit(numIgnores0,length(A0Trials));
[ignoresNAPct, ignoresNACI] = binofit(numIgnoresNA,length(AnaTrials));

ignoresByTrialType = 100*(horzcat(ignoresNAPct,ignores90Pct,ignores0Pct));
ignores90CI = ignores90Pct-ignores90CI;
ignores0CI = ignores0Pct-ignores0CI;
ignoresNACI = ignoresNAPct-ignoresNACI;
ignoresCI = vertcat(ignoresNACI,ignores90CI,ignores0CI);
% x = ["90deg Adaptor","0deg Adaptor","No Adaptor"];

%% Plot ignore rates for single mouse

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});

figure
plot(conditions,ignoresByTrialType, '-o')
hold on
title({'Ignore Rates per Adaptation Condition'; 'Subject: i1400'})
xlabel('Adaptation Condition')
ylabel('Ignore Rate (%)')
ylim([0,10])
hold off
%% Ignore Rates continued: Save multiple mice's outputs & plot

% ignoresByTrialType_1400 = ignoresByTrialType;
% ignoresByTrialType_601 = ignoresByTrialType;
ignoresByTrialType_1401 = ignoresByTrialType;

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});

figure
plot(conditions,ignoresByTrialType_1400, '-o')
hold on
plot(conditions,ignoresByTrialType_601, '-o')
plot(conditions,ignoresByTrialType_1401,'-o')
legend('i1400','i601','1401')
title({'Ignore Rates per Adaptation Condition'; 'Subjects: i1400,i601,i1401, 100ms Adapt'})
xlabel('Adaptation Condition')
ylabel('Ignore Rate (%)')
ylim([0,10])
hold off




%% Lapse Rate

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

LapseRateByTrialType = 100*(horzcat(AnaLapseRate,A90LapseRate,A0LapseRate));

%% Plot lapse rate for single animal

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});


figure
plot(conditions,LapseRateByTrialType, '-o')
hold on
title({'Lapse Rates per Adaptation Condition'; 'Subjects: i1401,i1402'})
xlabel('Adaptation Condition')
ylabel('Lapse Rate (%)')
ylim([0,10])
hold off

%% Lapse rates cont: Save for each animal & plot
% Lapse_Rates_601 = LapseRateByTrialType;
% Lapse_Rates_1400 = LapseRateByTrialType;
Lapse_Rates_1401 = LapseRateByTrialType;

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});


figure
plot(conditions,Lapse_Rates_1400, '-o')
hold on
plot(conditions,Lapse_Rates_601, '-o')
plot(conditions,Lapse_Rates_1401,'-o')
legend('i1400','i601','i1401')
title({'Lapse Rates per Adaptation Condition'; 'Subjects: i1400,i601,i1401, 100ms Adapt'})
xlabel('Adaptation Condition')
ylabel('Lapse Rate (%)')
ylim([0,10])
hold off

%% Delta change in ignore rate 

delta_ignores_1401 = [ignoresByTrialType_1401(2)/ignoresByTrialType_1401(1) ignoresByTrialType_1401(3)/ignoresByTrialType_1401(1)];
% delta_ignores_1402 = [ignoresByTrialType_1402(2)/ignoresByTrialType_1402(1) ignoresByTrialType_1402(3)/ignoresByTrialType_1402(1)];

conditions_90_0 = categorical(["90deg Adapter","0deg Adapter"]);
conditions_90_0 = reordercats(conditions_90_0,{'90deg Adapter','0deg Adapter'});

figure
plot(conditions_90_0,delta_ignores_1401,'-o')
hold on
plot(conditions_90_0,delta_ignores_1402,'-o')
legend('i1401','i1402')
title({'Delta Change in Ignore Rate in Presence of Adapter'})
xlabel('Adaptation Condition')
ylabel('Delta Ignore Rate (%)')
hold off 

%% Delta Change in lapse rate

delta_lapse_1401 = [Lapse_Rates_1401(2)/Lapse_Rates_1401(1) Lapse_Rates_1401(3)/Lapse_Rates_1401(1)];
% delta_lapse_1402 = [Lapse_Rates_1402(2)/Lapse_Rates_1402(1) Lapse_Rates_1402(3)/Lapse_Rates_1402(1)];

conditions_90_0 = categorical(["90deg Adapter","0deg Adapter"]);
conditions_90_0 = reordercats(conditions_90_0,{'90deg Adapter','0deg Adapter'});

figure
plot(conditions_90_0,delta_lapse_1401,'-o')
hold on
plot(conditions_90_0,delta_lapse_1402,'-o')
legend('i1401','i1402')
title({'Delta Change in Lapse Rates in Presence of Adapter'})
xlabel('Adaptation Condition')
ylabel('Delta Lapse Rate (%)')
hold off 










