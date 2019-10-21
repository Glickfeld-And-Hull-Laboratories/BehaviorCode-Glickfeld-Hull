% Avg reaction times per orientation, across adapt conditions 

%Goal here to combine data from different days into single matrices. 

successes = [];
incorrects = [];
lefttrials = [];
adaptori = [];
ori = [];
contrast = [];
DecisionTimeMs = [];

for i = 1:length(days_to_include)
    successes = [successes s(i).SIx];
    incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
    lefttrials = [lefttrials s(i).tLeftTrial];
    adaptori = [adaptori s(i).aGratingDirectionDeg];
    ori = [ori s(i).tGratingDirectionStart];
    contrast = [contrast s(i).aGratingContrast];
    DecisionTimeMs = [DecisionTimeMs s(i).tDecisionTimeMs];
end

DecisionTimeS = DecisionTimeMs./1000;

%Make matrix of ignores 
ignores = [];
for i = 1:length(successes)
    if successes(i)==0 & incorrects(i)==0
        ignores(i)=1;
    elseif successes(i)==1 | incorrects(i)==1
        ignores(i)=0;
    end
end
totaltrialsNI = length(successes)-(sum(ignores));

%NoAdapt indicates trials in the no-adaptor condition, where DoAdapt was
%turned on, but the contrast was set to 0
NoAdapt=nan(1,length(successes));
for itrial=1:length(successes)
    if contrast(itrial)==1
        NoAdapt(itrial)=0;
    elseif contrast(itrial)==0
        NoAdapt(itrial)=1;
    end
end

uniqueori = unique(ori);

%% Get rxn time per orientation - NA condition

AnaSuccesses = successes(NoAdapt==1);
AnaIncorrects = incorrects(NoAdapt==1);
AnaOriS = ori(NoAdapt==1 & successes==1);
AnaOriI = ori(NoAdapt==1 & incorrects==1);
AnaIgnores = ignores(NoAdapt==1);
AnaRxnTimeS = DecisionTimeS(NoAdapt==1 & successes==1);
AnaRxnTimeI = DecisionTimeS(NoAdapt==1 & incorrects==1);

AnaRxnTimeSPerOri=nan(1,length(uniqueori));
AnaRxnTimeIPerOri=nan(1,length(uniqueori));
AnaTrialSPerOri = nan(1,length(uniqueori));
AnaTrialIPerOri = nan(1,length(uniqueori));
AnaSemSPerOri = nan(1,length(uniqueori));
AnaSemIPerOri = nan(1,length(uniqueori));

for Anori=1:length(uniqueori)
    AnaRxnTimeSPerOri(Anori) = mean(AnaRxnTimeS(AnaOriS==uniqueori(Anori)));
    AnaRxnTimeIPerOri(Anori) = mean(AnaRxnTimeI(AnaOriI==uniqueori(Anori)));
    AnaTrialSPerOri(Anori) = sum(AnaSuccesses(AnaOriS==uniqueori(Anori)));
    AnaTrialIPerOri(Anori) = sum(AnaIncorrects(AnaOriI==uniqueori(Anori)));
    AnaSemSPerOri(Anori) = std(AnaRxnTimeS(AnaOriS==uniqueori(Anori)))/(sqrt(AnaTrialSPerOri(Anori)));
    AnaSemIPerOri(Anori) = std(AnaRxnTimeI(AnaOriI==uniqueori(Anori)))/(sqrt(AnaTrialIPerOri(Anori)));
end

% idx = logical([1 1 1 1 0 0 1 1 1 1]); 
% figure; errorbar(uniqueori(idx),AnaRxnTimeSPerOri(idx),AnaSemSPerOri(1,idx),AnaSemSPerOri(1,idx),'r','LineWidth',1.5)
% hold on
% errorbar(uniqueori(idx),AnaRxnTimeIPerOri(idx),AnaSemIPerOri(1,idx),AnaSemIPerOri(1,idx),'r','LineWidth',1.5)

%% get rxn time per orientation - 90deg condition

A90Successes = successes(adaptori==90 & NoAdapt==0);
A90Incorrects = incorrects(adaptori==90 & NoAdapt==0);
A90OriS = ori(adaptori==90 & NoAdapt==0 & successes==1);
A90OriI = ori(adaptori==90 & NoAdapt==0 & incorrects==1);
A90Ignores = ignores(adaptori==90 & NoAdapt==0);
A90RxnTimeS = DecisionTimeS(NoAdapt==0 & adaptori==90 & successes==1);
A90RxnTimeI = DecisionTimeS(NoAdapt==0 & adaptori==90 & incorrects==1);

A90RxnTimeSPerOri=nan(1,length(uniqueori));
A90RxnTimeIPerOri=nan(1,length(uniqueori));
A90TrialSPerOri = nan(1,length(uniqueori));
A90TrialIPerOri = nan(1,length(uniqueori));
A90SemSPerOri = nan(1,length(uniqueori));
A90SemIPerOri = nan(1,length(uniqueori));

for Anori=1:length(uniqueori)
    A90RxnTimeSPerOri(Anori) = mean(A90RxnTimeS(A90OriS==uniqueori(Anori)));
    A90RxnTimeIPerOri(Anori) = mean(A90RxnTimeI(A90OriI==uniqueori(Anori)));
    A90TrialSPerOri(Anori) = sum(A90Successes(A90OriS==uniqueori(Anori)));
    A90TrialIPerOri(Anori) = sum(A90Incorrects(A90OriI==uniqueori(Anori)));
    A90SemSPerOri(Anori) = std(A90RxnTimeS(A90OriS==uniqueori(Anori)))/(sqrt(A90TrialSPerOri(Anori)));
    A90SemIPerOri(Anori) = std(A90RxnTimeI(A90OriI==uniqueori(Anori)))/(sqrt(A90TrialIPerOri(Anori)));
end

% idx = logical([1 1 1 1 0 0 1 1 1 1]); 
% figure; errorbar(uniqueori(idx),A90RxnTimeSPerOri(idx),A90SemSPerOri(1,idx),A90SemSPerOri(1,idx),'r','LineWidth',1.5)
% hold on
% errorbar(uniqueori(idx),A90RxnTimeIPerOri(idx),A90SemIPerOri(1,idx),A90SemIPerOri(1,idx),'r','LineWidth',1.5)


%% Get rxn time per orientation - 0deg condition

A0Successes = successes(adaptori==0 & NoAdapt==0);
A0Incorrects = incorrects(adaptori==0 & NoAdapt==0);
A0OriS = ori(adaptori==0 & NoAdapt==0 & successes==1);
A0OriI = ori(adaptori==0 & NoAdapt==0 & incorrects==1);
A0Ignores = ignores(adaptori==0 & NoAdapt==0);
A0RxnTimeS = DecisionTimeS(NoAdapt==0 & adaptori==0 & successes==1);
A0RxnTimeI = DecisionTimeS(NoAdapt==0 & adaptori==0 & incorrects==1);

A0RxnTimeSPerOri=nan(1,length(uniqueori));
A0RxnTimeIPerOri=nan(1,length(uniqueori));
A0TrialSPerOri = nan(1,length(uniqueori));
A0TrialIPerOri = nan(1,length(uniqueori));
A0SemSPerOri = nan(1,length(uniqueori));
A0SemIPerOri = nan(1,length(uniqueori));

for Anori=1:length(uniqueori)
    A0RxnTimeSPerOri(Anori) = mean(A0RxnTimeS(A0OriS==uniqueori(Anori)));
    A0RxnTimeIPerOri(Anori) = mean(A0RxnTimeI(A0OriI==uniqueori(Anori)));
    A0TrialSPerOri(Anori) = sum(A0Successes(A0OriS==uniqueori(Anori)));
    A0TrialIPerOri(Anori) = sum(A0Incorrects(A0OriI==uniqueori(Anori)));
    A0SemSPerOri(Anori) = std(A0RxnTimeS(A0OriS==uniqueori(Anori)))/(sqrt(A0TrialSPerOri(Anori)));
    A0SemIPerOri(Anori) = std(A0RxnTimeI(A0OriI==uniqueori(Anori)))/(sqrt(A0TrialIPerOri(Anori)));
end

% idx = logical([1 1 1 1 0 0 1 1 1 1]); 
% figure; errorbar(uniqueori(idx),A0RxnTimeSPerOri(idx),A0SemSPerOri(1,idx),A0SemSPerOri(1,idx),'r','LineWidth',1.5)
% hold on
% errorbar(uniqueori(idx),A0RxnTimeIPerOri(idx),A0SemIPerOri(1,idx),A0SemIPerOri(1,idx),'r','LineWidth',1.5)


%%
% idx = logical([1 1 1 1 0 0 1 1 1 1]); 
% figure; errorbar(uniqueori(idx),AnaRxnTimeSPerOri(idx),AnaSemSPerOri(1,idx),AnaSemSPerOri(1,idx),'b','LineWidth',1.5)
% hold on
% errorbar(uniqueori(idx),A90RxnTimeSPerOri(idx),A90SemSPerOri(1,idx),A90SemSPerOri(1,idx),'r','LineWidth',1.5)
% errorbar(uniqueori(idx),A0RxnTimeSPerOri(idx),A0SemSPerOri(1,idx),A0SemSPerOri(1,idx),'g','LineWidth',1.5)
% title({'Mean Reaction Time per Orientation - Success Trials'; 'Paradigm: Flashing Adapt'; 'Subject: i1401'; 'n=3364 trials'})
% xlabel('Target Orientation (degrees)')
% ylabel('Mean Reaction Time (s)')
% legend('no adapter, n=949','90 deg adapter, n=1205','0 deg adapter, n=1210')
% hold off
% 
% 
% figure; errorbar(uniqueori(idx),AnaRxnTimeIPerOri(idx),AnaSemIPerOri(1,idx),AnaSemIPerOri(1,idx),'b','LineWidth',1.5)
% hold on
% errorbar(uniqueori(idx),A90RxnTimeIPerOri(idx),A90SemIPerOri(1,idx),A90SemIPerOri(1,idx),'r','LineWidth',1.5)
% errorbar(uniqueori(idx),A0RxnTimeIPerOri(idx),A0SemIPerOri(1,idx),A0SemIPerOri(1,idx),'g','LineWidth',1.5)
% title({'Mean Reaction Time per Orientation - Incorrect Trials'; 'Paradigm: Flashing Adapt'; 'Subject: i1401'; 'n=726 trials'})
% xlabel('Target Orientation (degrees)')
% ylabel('Mean Reaction Time (s)')
% legend('no adapter, n=188','90 deg adapter, n=252','0 deg adapter, n=286')
% hold off

figure; errorbar(uniqueori,AnaRxnTimeSPerOri,AnaSemSPerOri(1,:),AnaSemSPerOri(1,:),'b','LineWidth',1.5)
hold on
errorbar(uniqueori,A90RxnTimeSPerOri,A90SemSPerOri(1,:),A90SemSPerOri(1,:),'r','LineWidth',1.5)
errorbar(uniqueori,A0RxnTimeSPerOri,A0SemSPerOri(1,:),A0SemSPerOri(1,:),'g','LineWidth',1.5)
title({'Mean Reaction Time per Orientation - Success Trials'; 'Subject: i1402'; 'n=7089 trials'})
xlabel('Target Orientation (degrees)')
ylabel('Mean Reaction Time (s)')
legend('no adapter, n=4260','90 deg adapter, n=1449','0 deg adapter, n=1380')
hold off


figure; errorbar(uniqueori,AnaRxnTimeIPerOri,AnaSemIPerOri(1,:),AnaSemIPerOri(1,:),'b','LineWidth',1.5)
hold on
errorbar(uniqueori,A90RxnTimeIPerOri,A90SemIPerOri(1,:),A90SemIPerOri(1,:),'r','LineWidth',1.5)
errorbar(uniqueori,A0RxnTimeIPerOri,A0SemIPerOri(1,:),A0SemIPerOri(1,:),'g','LineWidth',1.5)
title({'Mean Reaction Time per Orientation - Incorrect Trials'; 'Subject: i1401'; 'n=1019 trials'})
xlabel('Target Orientation (degrees)')
ylabel('Mean Reaction Time (s)')
legend('no adapter, n=421','90 deg adapter, n=315','0 deg adapter, n=283')
hold off





