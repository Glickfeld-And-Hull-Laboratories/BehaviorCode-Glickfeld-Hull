%% Bootstrapping pSignifit to get avg slopes & CIs for eac adapt condition, each mouse
%Begin with days_to_include and s

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

%% Randomly sample 1000 trials from vars above using an index
trials = 1:length(successes);
idx = nan(100,1000);
slopies = nan(100,3);

for x = 1:100
    idx(x,:) = randsample(trials,1000);
end

for x = 1:100
    tLeftTrial_idx = tLeftTrial(idx(x,:));
    tLeftResponse_idx = tLeftResponse(idx(x,:));
    successes_idx = successes(idx(x,:));
    incorrects_idx = incorrects(idx(x,:));
    ori_idx = ori(idx(x,:));
    Aori_idx = Aori(idx(x,:));
    Acontrast_idx = Acontrast(idx(x,:));
    tIgnore_idx = tIgnore(idx(x,:));

    NoAdaptTrials = logical(~Acontrast_idx);
    Adapt90Trials = nan(1,length(Acontrast_idx));
    Adapt0Trials = nan(1,length(Acontrast_idx));

    for i=1:length(Acontrast_idx)
        if NoAdaptTrials(i)==0 & Aori_idx(i)==90
            Adapt90Trials(i)=1;
            Adapt0Trials(i)=0;
        elseif NoAdaptTrials(i)==0 & Aori_idx(i)==0
            Adapt90Trials(i)=0;
            Adapt0Trials(i)=1;
        elseif NoAdaptTrials(i)==1
            Adapt90Trials(i)=0;
            Adapt0Trials(i)=0;
        end
    end
    uniqueori = unique(ori_idx);

    NoAdaptTrialsLeft = tLeftResponse_idx(NoAdaptTrials==1 & tIgnore_idx==0);
    NoAdaptTrialsCorrect = successes_idx(NoAdaptTrials==1 & tIgnore_idx==0);
    NoAdaptTrialsIncorrect = incorrects_idx(NoAdaptTrials==1 & tIgnore_idx==0);
    NoAdaptOris = ori_idx(NoAdaptTrials==1 & tIgnore_idx==0);

    A90TrialsLeft = tLeftResponse_idx(Aori_idx==90 & NoAdaptTrials==0 & tIgnore_idx==0);
    A90TrialsCorrect = successes_idx(Aori_idx==90 & NoAdaptTrials==0 & tIgnore_idx==0);
    A90TrialsIncorrect = incorrects_idx(Aori_idx==90 & NoAdaptTrials==0 & tIgnore_idx==0);
    A90Oris = ori_idx(Aori_idx==90 & NoAdaptTrials==0 & tIgnore_idx==0);

    A0TrialsLeft = tLeftResponse_idx(Aori_idx==0 & NoAdaptTrials==0 & tIgnore_idx==0);
    A0TrialsCorrect = successes_idx(Aori_idx==0 & NoAdaptTrials==0 & tIgnore_idx==0);
    A0TrialsIncorrect = incorrects_idx(Aori_idx==0 & NoAdaptTrials==0 & tIgnore_idx==0);
    A0Oris = ori_idx(Aori_idx==0 & NoAdaptTrials==0 & tIgnore_idx==0);
    % For % LEFT! (Yes/No task) - array of # of left response, and total # trials for each ori
    NoAdaptLperOri=nan(1,length(uniqueori));
    A90LperOri=nan(1,length(uniqueori));
    A0LperOri=nan(1,length(uniqueori));
    
    NoAdaptNumperOri=nan(1,length(uniqueori));
    A90NumperOri=nan(1,length(uniqueori));
    A0NumperOri=nan(1,length(uniqueori));
    
    for o=1:length(uniqueori)
        NoAdaptLperOri(o) = sum(NoAdaptTrialsLeft(NoAdaptOris==uniqueori(o)));
        A90LperOri(o) = sum(A90TrialsLeft(A90Oris==uniqueori(o)));
        A0LperOri(o) = sum(A0TrialsLeft(A0Oris==uniqueori(o)));
        
        NoAdaptNumperOri(o) = sum(NoAdaptTrialsCorrect(NoAdaptOris==uniqueori(o))) + sum(NoAdaptTrialsIncorrect(NoAdaptOris==uniqueori(o)));
        A90NumperOri(o) = sum(A90TrialsCorrect(A90Oris==uniqueori(o))) + sum(A90TrialsIncorrect(A90Oris==uniqueori(o)));
        A0NumperOri(o) = sum(A0TrialsCorrect(A0Oris==uniqueori(o))) + sum(A0TrialsIncorrect(A0Oris==uniqueori(o)));
    end
    
    % idx = logical([1 1 1 1 0 0 1 1 1 1]);
    
    % output = vertcat(uniqueori(idx),NoAdaptLperOri(idx),NoAdaptNumperOri(idx))';
    % output90 = vertcat(uniqueori(idx),A90LperOri(idx),A90NumperOri(idx))';
    % output0 = vertcat(uniqueori(idx),A0LperOri(idx),A0NumperOri(idx))';
    output = vertcat(uniqueori,NoAdaptLperOri,NoAdaptNumperOri)';
    output90 = vertcat(uniqueori,A90LperOri,A90NumperOri)';
    output0 = vertcat(uniqueori,A0LperOri,A0NumperOri)';
    
    options = struct;
    options.sigmoidName    = 'logistic';
    options.expType        = 'YesNo';
    
    result=psignifit_EGB(output,options);
    result90=psignifit_EGB(output90,options);
    result0=psignifit_EGB(output0,options);
    
    Stimlevel = -45:0.5:45;
    slope = getSlope(result, Stimlevel);
    [slope_max, max_idx]=max(slope);
    slope90 = getSlope(result90, Stimlevel);
    [slope_max90, max_idx90]=max(slope90);
    slope0 = getSlope(result0, Stimlevel);
    [slope_max0, max_idx0]=max(slope0);

    slopies(x,:) = [slope_max,slope_max90,slope_max0];
end


%% Calculate CI's

avg_slopies = [mean(slopies(1,:)), mean(slopies(2,:)), mean(slopies(3,:))]
sd_slopies = [std(slopies(1,:)), std(slopies(2,:)), std(slopies(3,:))]

slopies_NA = sort(slopies(:,1));
slopies_NA_5 = slopies_NA(5);
slopies_NA_95 = slopies_NA(95);
CI_NA = vertcat(slopies_NA_5,slopies_NA_95);

slopies_90 = sort(slopies(:,2));
slopies_90_5 = slopies_90(5)
slopies_90_95 = slopies_90(95);
CI_90 = vertcat(slopies_90_5,slopies_90_95);

slopies_0 = sort(slopies(:,3));
slopies_0_5 = slopies_0(5);
slopies_0_95 = slopies_0(95);
CI_0 = vertcat(slopies_0_5,slopies_0_95);


%  CI = mean(x)+- t * (s / square(n))
%  t = 1.984, looked up in t-distribution table 
% where s is the standard deviation and n the sample size (= 100).   

% CI_slopies_hi = nan(1,3);
% CI_slopies_lo = nan(1,3);
% for j = 1:3
%     CI_slopies_hi(j) = avg_slopies(j) + 1.984*(sd_slopies(j)/sqrt(100));
%     CI_slopies_lo(j) = avg_slopies(j) - 1.984*(sd_slopies(j)/sqrt(100));
% end
% 
% alt_CI_slopies = vertcat(CI_slopies_lo,CI_slopies_hi)
%% Save stuff!

% slopies_1401 = slopies;
% avg_slopies_1401 = avg_slopies;
% sd_slopies_1401 = sd_slopies;
% CI_slopies_1401 = horzcat(CI_NA,CI_90,CI_0)

slopies_1402 = slopies;
avg_slopies_1402 = avg_slopies;
sd_slopies_1402 = sd_slopies;
CI_slopies_1402 = horzcat(CI_NA,CI_90,CI_0)

%% Plot slopes!
conditions = categorical(["No Adapter","90deg Adapter","0deg Adapter"]);
conditions = reordercats(conditions,{'No Adapter','90deg Adapter','0deg Adapter'});

figure
plot(conditions, max_slopes_1401,'-o','Color',[0 0.4470 0.7410])
hold on
f = errorbar(conditions,max_slopes_1401,(max_slopes_1401-CI_slopies_1401(1,:)),(CI_slopies_1401(2,:)-max_slopes_1401),'-o')
f.Color = [0 0.4470 0.7410];
plot(conditions,max_slopes_1402,'-o','Color',[0.8500 0.3250 0.0980])
e= errorbar(conditions, max_slopes_1402,(max_slopes_1402-CI_slopies_1402(1,:)),(CI_slopies_1402(2,:)-max_slopes_1402),'-o')
e.Color = [0.8500 0.3250 0.0980];
legend('i1401','i1401','i1402','i1402')
title({'Max Slope based on Psychometric Fits Across Conditions';'Subjects: i1401,i1402'})
xlabel('Adaptation Condition')
ylabel('Max Slope of Psychometric Curve')
ylim([0.015,.05])
hold off


    