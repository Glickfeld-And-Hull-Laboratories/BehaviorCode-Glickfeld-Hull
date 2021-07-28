% plot ignore rate by trial type

%%
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
    TGC = s_cropped_include(i).dGratingDiameter;
    dGratingDiameter = horzcat(dGratingDiameter, TGC); 
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

controls=find(tBlock2==0);
leds=find(tBlock2==1);
control_idx=zeros(1, length(tGratingContrast));
control_idx(controls)=1;
led_idx=tBlock2;
%% change contrast to diameter in this section for size analysis
target_ratio=round((tGratingContrast./dGratingContrast),3);
unq_target_ratios=flip(unique(target_ratio));
ratio1=unq_target_ratios(1); ratio2=unq_target_ratios(2); ratio3=unq_target_ratios(3); ratio4=unq_target_ratios(4);
ignore_trials=find(tIgnore);
left_ignore_idx=zeros(1,length(ignore_trials));
right_ignore_idx=zeros(1,length(ignore_trials));
for i=1:length(ignore_trials)
    islefttrial=tLeftTrial(ignore_trials(i));
    if islefttrial==1
        left_ignore_idx(i)=1;
    else
        right_ignore_idx(i)=1;
    end
    clear islefttrial
end


for i=1:length(ignore_trials) 
    targetratio(i)=tGratingContrast(ignore_trials(i))./dGratingContrast(ignore_trials(i));
end
targetratio=round(targetratio,3);
numtrials=length(tLeftTrial);

c1_left_ignore=length(intersect(find(left_ignore_idx), find(targetratio==ratio1)))/numtrials;
c2_left_ignore=length(intersect(find(left_ignore_idx), find(targetratio==ratio2)))/numtrials;
c3_left_ignore=length(intersect(find(left_ignore_idx), find(targetratio==ratio3)))/numtrials;
c4_left_ignore=length(intersect(find(left_ignore_idx), find(targetratio==ratio4)))/numtrials;

c1_right_ignore=length(intersect(find(right_ignore_idx), find(targetratio==ratio1)))/numtrials;
c2_right_ignore=length(intersect(find(right_ignore_idx), find(targetratio==ratio2)))/numtrials;  
c3_right_ignore=length(intersect(find(right_ignore_idx), find(targetratio==ratio3)))/numtrials;
c4_right_ignore=length(intersect(find(right_ignore_idx), find(targetratio==ratio4)))/numtrials;

plot_target_ratios=unq_target_ratios;
for i=1:4
   plot_target_ratios(i+4)=1/plot_target_ratios(i);
end
plot_target_ratios=sort(plot_target_ratios);
%unq_plot_ratios=unique(plot_target_ratios);
plot_mat(:,1)=plot_target_ratios';
plot_mat(:,2)=[c1_left_ignore, c2_left_ignore, c3_left_ignore, c4_left_ignore, c4_right_ignore, c3_right_ignore, c2_right_ignore, c1_right_ignore];
plot(plot_mat(:,1), plot_mat(:,2),'-blueo','LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor','blue', 'MarkerSize',5)
xlabel('Contrast Ratio'); ylabel('Percent of ignored trials')
title('i414 ignore rate by contrast ratio')
ax = gca

ax.XScale = 'log' 
xlim([0 102]) 
axis('square')

    
    