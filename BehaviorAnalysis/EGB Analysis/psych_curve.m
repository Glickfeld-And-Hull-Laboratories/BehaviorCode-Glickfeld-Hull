%% orientation psych curve creation

successes=s(1).SIx; %1 refers to day1, . is syntax to go in this column, SIx = success trials
incorrects=s(1).FIx; %FIx is incorrect trials only (not ignores)
lefttrials=s(1).tLeftTrial;
AdaptOri=s(1).aGratingDirectionDeg;
Ori=s(1).tGratingDirectionStart;


left_correct=intersect(find(successes), find(lefttrials));
right_correct=intersect(find(successes), find(~lefttrials)); %~ is not
left_incorrect=intersect(find(incorrects), find(lefttrials));

num_left_correct=length(left_correct);
num_right_correct=length(right_correct);
num_left_incorrect = length(left_incorrect);

Ltrials_only=find(lefttrials);
total_num_Ltrials=length(Ltrials_only); %successes serves as index of all trial bc 1 =success not success (incorrect and misses)
pct_left_correct=num_left_correct/total_num_Ltrials; %pct is percent
pct_right_correct=num_right_correct/total_num_Ltrials;
pct_left_incorrect=num_left_incorrect/total_num_Ltrials;

left_inc_corr_mat1=[pct_left_incorrect, pct_left_correct];

%figure(1)
%scatter(ori_deg,left_inc_corr_mat1,'filled')

   
%% 
successes=s(2).SIx; %1 refers to day1, . is syntax to go in this column, SIx = success trials
incorrects=s(2).FIx; %FIx is incorrect trials only (not ignores)
direction=s(2).tGratingDirectionStart;
lefttrials=s(2).tLeftTrial;

left_correct=intersect(find(successes), find(lefttrials));
right_correct=intersect(find(successes), find(~lefttrials)); %~ is not
left_incorrect=intersect(find(incorrects), find(lefttrials));

num_left_correct=length(left_correct)
num_right_correct=length(right_correct)
num_left_incorrect = length(left_incorrect)

Ltrials_only=find(lefttrials);
total_num_Ltrials=length(Ltrials_only); %successes serves as index of all trial bc 1 =success not success (incorrect and misses)
pct_left_correct=num_left_correct/total_num_Ltrials %pct is percent
pct_right_correct=num_right_correct/total_num_Ltrials
pct_left_incorrect=num_left_incorrect/total_num_Ltrials

left_inc_corr_mat2=[pct_left_incorrect, pct_left_correct]
%%
successes=s(3).SIx; %1 refers to day1, . is syntax to go in this column, SIx = success trials
incorrects=s(3).FIx; %FIx is incorrect trials only (not ignores)
direction=s(3).tGratingDirectionStart;
lefttrials=s(3).tLeftTrial;

left_correct=intersect(find(successes), find(lefttrials));
right_correct=intersect(find(successes), find(~lefttrials)); %~ is not
left_incorrect=intersect(find(incorrects), find(lefttrials));

num_left_correct=length(left_correct)
num_right_correct=length(right_correct)
num_left_incorrect = length(left_incorrect)

Ltrials_only=find(lefttrials);
total_num_Ltrials=length(Ltrials_only); %successes serves as index of all trial bc 1 =success not success (incorrect and misses)
pct_left_correct=num_left_correct/total_num_Ltrials %pct is percent
pct_right_correct=num_right_correct/total_num_Ltrials
pct_left_incorrect=num_left_incorrect/total_num_Ltrials

left_inc_corr_mat3=[pct_left_incorrect, pct_left_correct]

%%
%PUT PLOT CODE HERE
figure(2)
hold on
p1=plot(ori_deg,left_inc_corr_mat1,'-o')
p2=plot(ori_deg,left_inc_corr_mat2,'-o')
p3=plot(ori_deg,left_inc_corr_mat3,'-o')
p4=plot(ori_deg,left_inc_corr_mat4,'-o')
p5=plot(ori_deg,left_inc_corr_mat14,'-o')

%p1=plot(ori_deg,left_inc_corr_mat1,'-o')
%p2=plot(ori_deg,left_inc_corr_mat2,'-.r*')
%p3=plot(ori_deg,left_inc_corr_mat3,'--mo')
%p4=plot(ori_deg,left_inc_corr_mat4,':bs')
%p5=plot(ori_deg,left_inc_corr_mat5,'-.')

%p1.MarkerFaceColor='red'
%p1.MarkerEdgeColor='red'
