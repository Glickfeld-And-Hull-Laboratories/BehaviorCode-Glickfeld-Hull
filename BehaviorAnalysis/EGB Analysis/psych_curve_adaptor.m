%Run find_mult_files_for EGB_edited first to load all days
%Next run load_contrast_data_to_struct_forEGB_edited_602 to load all data
%to single structure 

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

%now have matrices of vars across all days, separated by adapt days and non-adapt days 

posoris = ori>0;
trialtypetest = sum(lefttrials-posoris);
Aposoris = Aori>0;
Atrialtypetest = sum(Alefttrials-Aposoris);
%+oris = Ltrials, -oris = Rtrials. Test that this aligns w tlefttrial var
%trialtypetest should =0

%%
%psych curve, percent Left 

% uniqueori = unique(ori);
% success=nan(1,length(uniqueori));
% pctL=[];
% 
% for nori=1:length(uniqueori)
%     success(nori)= sum(successes(ori==uniqueori(nori)));
%     alloritrials(nori) = length(successes(ori==uniqueori(nori)));
%     if uniqueori(nori)>0
%         pctL(nori) = success(nori)./alloritrials(nori);
%     elseif uniqueori(nori)<0
%         pctL(nori) = 1-(success(nori)./alloritrials(nori));
%     end
% end

%need to do same for adapt condition
%Separate out by contrast for interleaved adapt/no adapt


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
    A90success(Anori) = sum(A90Trials(A90Oris==Auniqueori(Anori)));
    A90alloritrials(Anori) = length(A90Trials(A90Oris==Auniqueori(Anori)));
    if Auniqueori(Anori)>0
        A90pctL(Anori) = A90success(Anori)./A90alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A90pctL(Anori) = 1-(A90success(Anori)./A90alloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    A0success(Anori) = sum(A0Trials(A0Oris==Auniqueori(Anori)));
    A0alloritrials(Anori) = length(A0Trials(A0Oris==Auniqueori(Anori)));
    if Auniqueori(Anori)>0
        A0pctL(Anori) = A0success(Anori)./A0alloritrials(Anori);
    elseif Auniqueori(Anori)<0
        A0pctL(Anori) = 1-(A0success(Anori)./A0alloritrials(Anori));
    end
end

for Anori=1:length(Auniqueori)
    Anasuccess(Anori) = sum(AnaTrials(AnaOris==Auniqueori(Anori)));
    Anaalloritrials(Anori) = length(AnaTrials(AnaOris==Auniqueori(Anori)));
    if Auniqueori(Anori)>0
        AnapctL(Anori) = Anasuccess(Anori)./Anaalloritrials(Anori);
    elseif Auniqueori(Anori)<0
        AnapctL(Anori) = 1-(Anasuccess(Anori)./Anaalloritrials(Anori));
    end
end

% figure;plot(uniqueori, pctL,'k')
% hold on
figure;plot(Auniqueori, A90pctL,'r')
hold on
plot(Auniqueori, A0pctL,'g')
plot(Auniqueori, AnapctL,'b')
%%
%pct correct 

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
    A90success(Anori) = sum(A90Trials(A90Oris==Auniqueori(Anori)));
    A90alloritrials(Anori) = length(A90Trials(A90Oris==Auniqueori(Anori)));
    A90pct(Anori) = A90success(Anori)./A90alloritrials(Anori);
end

for Anori=1:length(Auniqueori)
    A0success(Anori) = sum(A0Trials(A0Oris==Auniqueori(Anori)));
    A0alloritrials(Anori) = length(A0Trials(A0Oris==Auniqueori(Anori)));
    A0pct(Anori) = A0success(Anori)./A0alloritrials(Anori);
end

for Anori=1:length(Auniqueori)
    Anasuccess(Anori) = sum(AnaTrials(AnaOris==Auniqueori(Anori)));
    Anaalloritrials(Anori) = length(AnaTrials(AnaOris==Auniqueori(Anori)));
    Anapct(Anori) = Anasuccess(Anori)./Anaalloritrials(Anori);
end

% figure;plot(uniqueori, pct,'k')
% hold on
figure;plot(Auniqueori, A90pct,'r')
hold on
plot(Auniqueori, A0pct,'g')
plot(Auniqueori, Anapct,'b')



