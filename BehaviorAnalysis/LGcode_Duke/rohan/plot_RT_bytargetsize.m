% plot reaction time by target size or contrast

SIx = s_cropped_include.SIx;
for i = 2:length(s_cropped_include); 
    success = s_cropped_include(i).SIx; 
    SIx = horzcat(SIx, success);
end 

FIx = s_cropped_include.FIx;
for i = 2:length(s_cropped_include); 
    failure = s_cropped_include(i).FIx; 
    FIx = horzcat(FIx, failure); 
end 

tLeftTrial = s_cropped_include.tLeftTrial;
for i = 2:length(s_cropped_include);
    TLT = s_cropped_include(i).tLeftTrial;
    tLeftTrial = horzcat(tLeftTrial, TLT); 
end

tGratingContrast = s_cropped_include.tGratingContrast;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).tGratingContrast;
    tGratingContrast = horzcat(tGratingContrast, TGC); 
end 

dGratingContrast = s_cropped_include.dGratingContrast;
for i = 2:length(s_cropped_include); 
    DGC = s_cropped_include(i).dGratingContrast;
    dGratingContrast = horzcat(dGratingContrast, DGC); 
end 

tGratingDiameter = s_cropped_include.tGratingDiameter;
for i = 2:length(s_cropped_include); 
    TGC = s_cropped_include(i).tGratingDiameter;
    tGratingDiameter = horzcat(tGratingDiameter, TGC); 
end 

dGratingDiameter = s_cropped_include.dGratingDiameter;
for i = 2:length(s_cropped_include); 
    DGC = s_cropped_include(i).dGratingDiameter;
    dGratingDiameter = horzcat(dGratingDiameter, DGC); 
end 

tRightResponse = s_cropped_include.tRR;
for i = 2:length(s_cropped_include); 
    tR_resp = s_cropped_include(i).tRR;
    tRightResponse = horzcat(tRightResponse, tR_resp); 
end 

tIgnore = s_cropped_include.tIgnore;
for i = 2:length(s_cropped_include); 
    Ignore = s_cropped_include(i).tIgnore;
    tIgnore = horzcat(tIgnore, Ignore); 
end 

tBlock2 = s_cropped_include.tBlock2; 
for i = 2:length(s_cropped_include); 
    TBT = s_cropped_include(i).tBlock2; 
    tBlock2 = horzcat(tBlock2, TBT); 
end 

decisiontimes= s_cropped_include.tDecisionTimeMs
for i= 2:length(s_cropped_include);
    decisiontimesval=s_cropped_include(i).tDecisionTimeMs;
    decisiontimes=horzcat(decisiontimes, decisiontimesval);
end


controls=find(tBlock2==0);
leds=find(tBlock2==1);
control_idx=zeros(1, length(tGratingContrast));
control_idx(controls)=1;
led_idx=tBlock2;

%%
unqtargets=unique(tGratingDiameter);
unqtargets=round(unqtargets,4);
t1=unqtargets(1); t2=unqtargets(2); t3=unqtargets(3); t4=unqtargets(4);

resp_trials=find(~tIgnore);
resp_idx=zeros(1,length(tGratingContrast));
resp_idx(resp_trials)=1;

tGratingDiameter=round(tGratingDiameter, 4);

t1_trials=intersect(find(resp_idx), find(tGratingDiameter==t1));
t2_trials=intersect(find(resp_idx), find(tGratingDiameter==t2));
t3_trials=intersect(find(resp_idx), find(tGratingDiameter==t3));
t4_trials=intersect(find(resp_idx), find(tGratingDiameter==t4));

t1_dt = decisiontimes(t1_trials);
t2_dt= decisiontimes(t2_trials);
t3_dt=decisiontimes(t3_trials);
t4_dt=decisiontimes(t4_trials);
t1_mean_dt=mean(t1_dt);
t2_mean_dt=mean(t2_dt);
t3_mean_dt=mean(t3_dt);
t4_mean_dt=mean(t4_dt);
mean_dectimes=[t1_mean_dt, t2_mean_dt, t3_mean_dt, t4_mean_dt];

figure;
hold on
scatter(unqtargets, mean_dectimes,'filled')
xlabel('Target Stimulus Diameter (Deg)'); ylabel('Reaction Time (ms)')
title('i433 mean reaction time by target stimulus size')




