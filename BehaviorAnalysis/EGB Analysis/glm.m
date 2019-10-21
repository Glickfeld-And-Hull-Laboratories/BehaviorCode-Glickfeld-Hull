%load input of vars & data, single day or collapsed 
%Which vars am I interested in eval the weights contributing to var? 
%trialOutcomeCell
%tLeftTrial
%Control: tGratingContrast - same across trials 
%tDecisionTimeMs
%tNTrialsCompleted
%tConsecErrors
%tConsecCorrects
%aGratingContrast
%aGratingDirectionDeg


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
        trialOutcomeCell(i)=2;
    elseif ignores(i)==0
        trialOutcomeCell(i)=trialOutcomeCell(i);
    end
end

SimpleTrialOutcomeCell = successes;

Ori = double(cell2mat(input.tGratingDirectionStart));
tLeftTrial = double(cell2mat(input.tLeftTrial));
TargetContrastControl = double(cell2mat(input.tGratingContrast));
RxnTime = double(cell2mat(input.tDecisionTimeMs));
NTrialsCompleted = double(cell2mat(input.tNTrialsCompleted));
ConsecErrors = double(cell2mat(input.tConsecErrors));
ConsecCorrects = double(cell2mat(input.tConsecCorrects));
aGratingContrast = double(cell2mat(input.aGratingContrast));
aGratingDirectionDeg = double(cell2mat(input.aGratingDirectionDeg));

All_Vars = vertcat(TargetContrastControl,Ori,tLeftTrial,RxnTime,...
    NTrialsCompleted,ConsecErrors,ConsecCorrects,aGratingContrast,aGratingDirectionDeg);
All_Vars = All_Vars';

Some_Vars = vertcat(TargetContrastControl,Ori);
Some_Vars = Some_Vars';


N_Trials = [1:length(successes)];
Y_Var = vertcat(SimpleTrialOutcomeCell,N_Trials)';

b = glmfit(All_Vars,Y_Var,'binomial')

