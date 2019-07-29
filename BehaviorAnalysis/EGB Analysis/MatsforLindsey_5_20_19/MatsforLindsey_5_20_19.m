% 2 matrices for each of the three mice on the orientation discrimination task with the 100 ms adapter: 
% 1) the number of left choices (so it increases if you have a -ori -> + ori x-axis) for each orientation in the three adaptation conditions 
% (so this should be nOri x 3) and 
% 2) the number of total trials for each orientation in the three adaptation conditions? Thanks!

%combine data from different days into single matrices. 
tLeftTrial = [];
tLeftResponse = [];
successes = [];
incorrects = [];
ori = [];
Aori = [];
Acontrast = [];

for i=1:length(days_to_include)
    tLeftTrial = [tLeftTrial s(i).tLeftTrial];
    tLeftResponse = [tLeftResponse s(i).tLeftResponse];
    successes = [successes s(i).SIx];
    incorrects = [incorrects s(i).FIx]; %FIx is incorrect trials only (not ignores)
    ori = [ori s(i).tGratingDirectionStart];
    Aori = [Aori s(i).aGratingDirectionDeg];
    Acontrast = [Acontrast s(i).aGratingContrast];
end

%1 mat for ea adapt condition, 1s=that adapt condition, 0s=not
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

%For loop for each orientation to output #L per ori. Do vars for each adapt
%condition
uniqueori = unique(ori);

NoAdaptTrialsLeft = tLeftResponse(NoAdaptTrials==1);
NoAdaptOris = ori(NoAdaptTrials==1);

A90TrialsLeft = tLeftResponse(Aori==90 & NoAdaptTrials==0);
A90Oris = ori(Aori==90 & NoAdaptTrials==0);

A0TrialsLeft = tLeftResponse(Aori==0 & NoAdaptTrials==0);
A0Oris = ori(Aori==0 & NoAdaptTrials==0);

NoAdaptLperOri=nan(1,length(uniqueori)); %array of # of total trials left response for each ori
A90LperOri=nan(1,length(uniqueori));
A0LperOri=nan(1,length(uniqueori));

for ori=1:length(uniqueori)
    NoAdaptLperOri(ori) = sum(NoAdaptTrialsLeft(NoAdaptOris==uniqueori(ori)));
    A90LperOri(ori) = sum(A90TrialsLeft(A90Oris==uniqueori(ori)));
    A0LperOri(ori) = sum(A0TrialsLeft(A0Oris==uniqueori(ori)));
end

NumLperOri_NA_0_90 = vertcat(NoAdaptLperOri,A0LperOri,A90LperOri);
save('NumLperOri601.mat','NumLperOri_NA_0_90')

    
    
    
    
    
    
    
