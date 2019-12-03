%Ignore rate across days 

%Run find_mult_files_for EGB_edited first to load all days
%Next run load_contrast_data_to_struct_forEGB_edited_602 to load all data
%to single structure 

%First need to separate out adapt and non-adapt days 
% doadapt = []; 
% for i = 1:length(days_to_include)
%     doadapt = [doadapt s(i).doAdapt];
% end

%Matrices of important metrics on each day
SuccessesPerDay = [];
IncorrectsPerDay = [];
TrialsPerDay = [];

for i=1:length(days_to_include)
    SuccessesPerDay = [SuccessesPerDay sum(s(i).SIx)];
    IncorrectsPerDay = [IncorrectsPerDay sum(s(i).FIx)];
    TrialsPerDay = [TrialsPerDay length(s(i).SIx)];
end

IgnoresPerDay = TrialsPerDay - (SuccessesPerDay + IncorrectsPerDay);
IgnoreRatePerDay = IgnoresPerDay./TrialsPerDay;

f = fit(days_to_include', IgnoreRatePerDay','exp1')
figure; plot(f, days_to_include', IgnoreRatePerDay')
hold on

title({'Ignore Rate Across Days in Early Training Stage'; 'Subject: i1400'; '11/21 - 11/30'})
xlabel('Day of Training')
ylabel('Ignore Rate')
xlim([1,7])
hold off 

%Have to edit by hand b/c of fit function 

%Make same plot but do ignore rate for every 100 trials 
successes = [];
incorrects = [];
lefttrials = [];
ori = [];

for i = 1:length(days_to_include)
        successes = [successes s(i).SIx];
        incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
        lefttrials = [lefttrials s(i).tLeftTrial];
        ori = [ori s(i).tGratingDirectionStart];
end

%Make matrix of ignores 
ignores = [];
for i = 1:length(successes)
    if successes(i)==0 & incorrects(i)==0
        ignores(i)=1;
    elseif successes(i)==1 | incorrects(i)==1
        ignores(i)=0;
    end
end

%Make matrix of ignore rate per 100 trials 
ignoresPer100 = [];
for i=1:(floor(length(ignores)./100))
    ignoresPer100(i) = sum(ignores((1+((i-1)*100)):(100*i)));
end
% N=(round(length(ignores)./100));    
% temp=[];
% temp=reshape(ignores(1:N*100),100,N);
% ignoresPer100_new= sum(temp,1);
% 

AllTrials = floor(length(ignores)./100)*100;
Each100trials = [];

for i = 1:AllTrials/100
    Each100trials(i) = 100*i;
end

f = fit(Each100trials', ignoresPer100','exp1')
figure; plot(f, Each100trials', ignoresPer100')
hold on

title({'Ignore Rate per 100 Trials in Early Training Stage'; 'Subject: i1400'; '11/21 - 11/30'})
xlabel('Group of 100 Trials')
ylabel('# of Ignores per 100 Trials')
xlim([0,2200])
hold off 




