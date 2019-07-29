%Desired output is input for psignifit:  "pass your data in the n x 3 matrix of the form:
%       [x-value, number correct, number of trials]"

%combine data from different days into single matrices. 
tLeftTrial = [];
tLeftResponse = [];
successes = [];
incorrects = [];
ori = [];
Aori = [];
Acontrast = [];
tIgnore = [];

for i=1:length(days_to_include)
    tLeftTrial = [tLeftTrial s(i).tLeftTrial];
    tLeftResponse = [tLeftResponse s(i).tLeftResponse];
    successes = [successes s(i).SIx];
    incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
    ori = [ori s(i).tGratingDirectionStart];
    Aori = [Aori s(i).aGratingDirectionDeg];
    Acontrast = [Acontrast s(i).aGratingContrast];
    tIgnore = [tIgnore s(i).tIgnore];
end



%%
%1 mat for ea adapt condition, 1;s=that adapt condition, 0;s=not
NoAdaptTrials = logical(~Acontrast);
Adapt90Trials = nan(1,length(Acontrast));
Adapt0Trials = nan(1,length(Acontrast));

for i=1:length(Acontrast)
    if NoAdaptTrials(i)==0 & Aori(i)==90
        Adapt90Trials(i)=1;
        Adapt0Trials(i)=0;
    elseif NoAdaptTrials(i)==0 & Aori(i)==0
        Adapt90Trials(i)=0;
        Adapt0Trials(i)=1;
    elseif NoAdaptTrials(i)==1
        Adapt90Trials(i)=0;
        Adapt0Trials(i)=0;
    end
end

%%
%For loop for each orientation to output NumCorrect and NumTrials (correct + incorrect) per ori. Do vars for each adapt
%condition
uniqueori = unique(ori);

NoAdaptTrialsLeft = tLeftResponse(NoAdaptTrials==1 & tIgnore==0);
NoAdaptTrialsCorrect = successes(NoAdaptTrials==1 & tIgnore==0);
NoAdaptTrialsIncorrect = incorrects(NoAdaptTrials==1 & tIgnore==0);
NoAdaptOris = ori(NoAdaptTrials==1 & tIgnore==0);

A90TrialsLeft = tLeftResponse(Aori==90 & NoAdaptTrials==0 & tIgnore==0);
A90TrialsCorrect = successes(Aori==90 & NoAdaptTrials==0 & tIgnore==0);
A90TrialsIncorrect = incorrects(Aori==90 & NoAdaptTrials==0 & tIgnore==0);
A90Oris = ori(Aori==90 & NoAdaptTrials==0 & tIgnore==0);

A0TrialsLeft = tLeftResponse(Aori==0 & NoAdaptTrials==0 & tIgnore==0);
A0TrialsCorrect = successes(Aori==0 & NoAdaptTrials==0 & tIgnore==0);
A0TrialsIncorrect = incorrects(Aori==0 & NoAdaptTrials==0 & tIgnore==0);
A0Oris = ori(Aori==0 & NoAdaptTrials==0 & tIgnore==0);

%%
%array of # of left response, and total # trials for each ori
NoAdaptLperOri=nan(1,length(uniqueori)); 
A90LperOri=nan(1,length(uniqueori));
A0LperOri=nan(1,length(uniqueori));

NoAdaptNumperOri=nan(1,length(uniqueori)); 
A90NumperOri=nan(1,length(uniqueori));
A0NumperOri=nan(1,length(uniqueori));

for ori=1:length(uniqueori) 
    
    
    NoAdaptLperOri(ori) = sum(NoAdaptTrialsLeft(NoAdaptOris==uniqueori(ori)));
    A90LperOri(ori) = sum(A90TrialsLeft(A90Oris==uniqueori(ori)));
    A0LperOri(ori) = sum(A0TrialsLeft(A0Oris==uniqueori(ori)));
    
    NoAdaptNumperOri(ori) = sum(NoAdaptTrialsCorrect(NoAdaptOris==uniqueori(ori))) + sum(NoAdaptTrialsIncorrect(NoAdaptOris==uniqueori(ori)));
    A90NumperOri(ori) = sum(A90TrialsCorrect(A90Oris==uniqueori(ori))) + sum(A90TrialsIncorrect(A90Oris==uniqueori(ori)));
    A0NumperOri(ori) = sum(A0TrialsCorrect(A0Oris==uniqueori(ori))) + sum(A0TrialsIncorrect(A0Oris==uniqueori(ori)));
end

idx = logical([1 1 1 1 0 0 1 1 1 1]);

output = vertcat(uniqueori(idx),NoAdaptLperOri(idx),NoAdaptNumperOri(idx))';
output90 = vertcat(uniqueori(idx),A90LperOri(idx),A90NumperOri(idx))';
output0 = vertcat(uniqueori(idx),A0LperOri(idx),A0NumperOri(idx))';
options = struct;
options.sigmoidName    = 'logistic';
options.expType        = 'YesNo';
%%
cd('C:\Users\emily2\Documents\Repositories\BehaviorCode-Glickfeld-Hull\BehaviorAnalysis\EGB Analysis\psignifit-master\psignifit-master')

result=psignifit(output,options);
result90=psignifit(output90,options);
result0=psignifit(output0,options);

plotOptions= struct;
plotOptions.dataColor      = [1,0,0];
plotOptions.lineColor      = [0.6350, 0.0780, 0.1840]; 

plotOptions2= struct;
plotOptions2.dataColor      = [0.4660, 0.6740, 0.1880];
plotOptions2.lineColor      = [0, 0.5, 0];


figure; plotPsych(result)
hold on
plotPsych(result90,plotOptions)
plotPsych2(result0,plotOptions2)
title({'Logistic Fit for Psychometric Curve'; 'Paradigm: Flashing Adapt'; 'Subject: i1401'})
xlabel('Target Orientation (degrees)')
ylabel('% Left Response')
legend('no adapter, n=1137','90 deg adapter, n=1457','0 deg adapter, n=1496')
hold off

Stimlevel= -45:0.5:45;
slope = getSlope(result, Stimlevel);
[slope_max, max_idx]=max(slope);
slope90 = getSlope(result90, Stimlevel);
[slope_max90, max_idx90]=max(slope90);
slope0 = getSlope(result0, Stimlevel);
[slope_max0, max_idx0]=max(slope0);

figure
plot(Stimlevel, slope,'b')
hold on
plot(Stimlevel, slope90, 'r')
plot(Stimlevel, slope0, 'g')
ylim([0,.05]);
title({'Slopes'; 'Paradigm: Flashing Adapt'; 'Subject: i1401'})
xlabel('Target Orientation (degrees)')
ylabel('Slope of Psych Curve Fit')
legend('no adapter, n=1137','90 deg adapter, n=1457','0 deg adapter, n=1496')
vline(Stimlevel(max_idx),'b:')
hline(slope_max,'b:')
vline(Stimlevel(max_idx90),'r:')
hline(slope_max90,'r:')
vline(Stimlevel(max_idx0),'g:')
hline(slope_max0,'g:')
hold off

cd('Z:\home\emily\2018-19\Figures\Joint_Lab_mtg_06.13.19\4s_Adapt');
save('4sec_Adapt_psignifit_workspace');

%%
%%Quantify difference in max slope
conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});
max_slopes = horzcat(slope_max,slope_max90,slope_max0);
figure
b = bar(conditions,max_slopes,'FaceColor',[1 1 1],'LineWidth',1.5)
hold on
b.EdgeColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
title({'Slopes of Psychometric Fits Across Conditions'; 'Paradigm: 4s Static Adapt'; 'Subjects: i601,i1401'; '12/03-02/07'})
ylabel('Slope of Psych Curve Fit')
hold off

conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});
max_slopes = horzcat(slope_max,slope_max90,slope_max0);
figure
b = bar(conditions,max_slopes,'FaceColor',[1 1 1],'LineWidth',1.5)
hold on
b.EdgeColor = 'flat';
b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];
b.CData(3,:) = [0 1 0];
title({'Slopes of Psychometric Fits Across Conditions'; 'Paradigm: 4s Static Adapt'; 'Subjects: i601,i1401'; '12/03-02/07'})
ylabel('Slope of Psych Curve Fit')
hold off

