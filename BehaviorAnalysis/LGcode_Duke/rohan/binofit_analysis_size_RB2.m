%% bootsrap mean distribution analysis

% get concatenated data
%isolate LED and control trials if not already isolated
%get random number generator to pull ~250 trials and get mean percent
%correct?

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

%% CONTRAST SORTING

contrastratio=tGratingDiameter./dGratingDiameter;
contrastratio=round(contrastratio,4);
contrast1_rightidx=zeros(1, length(tGratingContrast)); %declare contrast idx arrays
contrast1_leftidx=contrast1_rightidx;
contrast2_rightidx=contrast1_rightidx;
contrast2_leftidx=contrast1_rightidx;
contrast3_rightidx=contrast1_rightidx;
contrast3_leftidx=contrast1_rightidx;
contrast4_rightidx=contrast1_rightidx;
contrast4_leftidx=contrast1_rightidx;

unqcontrast=unique(round(contrastratio,4)); %CONTRAST SORTED BY SIDE
contrast1_righttrials=intersect(find(contrastratio==unqcontrast(4)), find(~tLeftTrial));
contrast1_rightidx(contrast1_righttrials)=1;
contrast1_lefttrials=intersect(find(contrastratio==unqcontrast(4)), find(tLeftTrial));
contrast1_leftidx(contrast1_lefttrials)=1;
contrast2_righttrials=intersect(find(contrastratio==unqcontrast(3)), find(~tLeftTrial));
contrast2_rightidx(contrast2_righttrials)=1;
contrast2_lefttrials=intersect(find(contrastratio==unqcontrast(3)), find(tLeftTrial));
contrast2_leftidx(contrast2_lefttrials)=1;
contrast3_righttrials=intersect(find(contrastratio==unqcontrast(2)), find(~tLeftTrial));
contrast3_rightidx(contrast3_righttrials)=1;
contrast3_lefttrials=intersect(find(contrastratio==unqcontrast(2)), find(tLeftTrial));
contrast3_leftidx(contrast3_lefttrials)=1;
contrast4_righttrials=intersect(find(contrastratio==unqcontrast(1)), find(~tLeftTrial));
contrast4_rightidx(contrast4_righttrials)=1;
contrast4_lefttrials=intersect(find(contrastratio==unqcontrast(1)), find(tLeftTrial));
contrast4_leftidx(contrast4_lefttrials)=1;
c1=unqcontrast(4);
c2=unqcontrast(3);
c3=unqcontrast(2);
c4=unqcontrast(1);

%% sort contrast sides into control and LED
c1_right_ctrl_idx=zeros(1,length(tGratingContrast));
c1_left_ctrl_idx=c1_right_ctrl_idx;
c2_right_ctrl_idx=c1_right_ctrl_idx;
c2_left_ctrl_idx=c1_right_ctrl_idx;
c3_right_ctrl_idx=c1_right_ctrl_idx;
c3_left_ctrl_idx=c1_right_ctrl_idx;
c4_right_ctrl_idx=c1_right_ctrl_idx;
c4_left_ctrl_idx=c1_right_ctrl_idx;

c1_right_ctrl=intersect(find(contrast1_rightidx), find(control_idx));
c1_right_ctrl_idx(c1_right_ctrl)=1;
c1_left_ctrl=intersect(find(contrast1_leftidx), find(control_idx));
c1_left_ctrl_idx(c1_left_ctrl)=1;
c2_right_ctrl=intersect(find(contrast2_rightidx), find(control_idx));
c2_right_ctrl_idx(c2_right_ctrl)=1;
c2_left_ctrl=intersect(find(contrast2_leftidx), find(control_idx));
c2_left_ctrl_idx(c2_left_ctrl)=1;
c3_right_ctrl=intersect(find(contrast3_rightidx), find(control_idx));
c3_right_ctrl_idx(c3_right_ctrl)=1;
c3_left_ctrl=intersect(find(contrast3_leftidx), find(control_idx));
c3_left_ctrl_idx(c3_left_ctrl)=1;
c4_right_ctrl=intersect(find(contrast4_rightidx), find(control_idx));
c4_right_ctrl_idx(c4_right_ctrl)=1;
c4_left_ctrl=intersect(find(contrast4_leftidx), find(control_idx));
c4_left_ctrl_idx(c4_left_ctrl)=1;

c1_right_ctrl_correct_idx=zeros(1,length(tGratingContrast));
c1_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c2_right_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c2_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c3_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c3_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c4_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;
c4_left_ctrl_correct_idx=c1_right_ctrl_correct_idx;

% find correct trials of controls already sorted by contrast and side
c1_right_ctrl_correct=intersect(find(c1_right_ctrl_idx), find(SIx));
c1_left_ctrl_correct=intersect(find(c1_left_ctrl_idx), find(SIx));
c1_right_ctrl_correct_idx(c1_right_ctrl_correct)=1;
c1_left_ctrl_correct_idx(c1_left_ctrl_correct)=1;
c2_right_ctrl_correct=intersect(find(c2_right_ctrl_idx), find(SIx));
c2_left_ctrl_correct=intersect(find(c2_left_ctrl_idx), find(SIx));
c2_right_ctrl_correct_idx(c2_right_ctrl_correct)=1;
c2_left_ctrl_correct_idx(c2_left_ctrl_correct)=1;
c3_right_ctrl_correct=intersect(find(c3_right_ctrl_idx), find(SIx));
c3_left_ctrl_correct=intersect(find(c3_left_ctrl_idx), find(SIx));
c3_right_ctrl_correct_idx(c3_right_ctrl_correct)=1;
c3_left_ctrl_correct_idx(c3_left_ctrl_correct)=1;
c4_right_ctrl_correct=intersect(find(c4_right_ctrl_idx), find(SIx));
c4_left_ctrl_correct=intersect(find(c4_left_ctrl_idx), find(SIx));
c4_right_ctrl_correct_idx(c4_right_ctrl_correct)=1;
c4_left_ctrl_correct_idx(c4_left_ctrl_correct)=1;
%num of trials is ex: sum(c4_left_ctrl_idx) so binofit(


%% repeat above process with LED trials
c1_right_led_idx=zeros(1,length(tGratingContrast));
c1_left_led_idx=c1_right_led_idx;
c2_right_led_idx=c1_right_led_idx;
c2_left_led_idx=c1_right_led_idx;
c3_right_led_idx=c1_right_led_idx;
c3_left_led_idx=c1_right_led_idx;
c4_right_led_idx=c1_right_led_idx;
c4_left_led_idx=c1_right_led_idx;
c1_right_led_correct_idx=c1_right_led_idx;
c1_left_led_correct_idx=c1_right_led_idx;
c2_right_led_correct_idx=c1_right_led_idx;
c2_left_led_correct_idx=c1_right_led_idx;
c3_right_led_correct_idx=c1_right_led_idx;
c3_left_led_correct_idx=c1_right_led_idx;
c4_right_led_correct_idx=c1_right_led_idx;
c4_left_led_correct_idx=c1_right_led_idx;

c1_right_led=intersect(find(contrast1_rightidx), find(led_idx));
c1_left_led=intersect(find(contrast1_leftidx), find(led_idx));
c1_right_led_idx(c1_right_led)=1;
c1_left_led_idx(c1_left_led)=1;
c2_right_led=intersect(find(contrast2_rightidx), find(led_idx));
c2_left_led=intersect(find(contrast2_leftidx), find(led_idx));
c2_right_led_idx(c2_right_led)=1;
c2_left_led_idx(c2_left_led)=1;
c3_right_led=intersect(find(contrast3_rightidx), find(led_idx));
c3_left_led=intersect(find(contrast3_leftidx), find(led_idx));
c3_right_led_idx(c3_right_led)=1;
c3_left_led_idx(c3_left_led)=1;
c4_right_led=intersect(find(contrast1_rightidx), find(led_idx));
c4_left_led=intersect(find(contrast4_leftidx), find(led_idx));
c4_right_led_idx(c4_right_led)=1;
c4_left_led_idx(c4_left_led)=1;

c1_right_led_correct=intersect(find(c1_right_led_idx), find(SIx));
c1_left_led_correct=intersect(find(c1_left_led_idx), find(SIx));
c1_right_led_correct_idx(c1_right_led_correct)=1;
c1_left_led_correct_idx(c1_left_led_correct)=1;
c2_right_led_correct=intersect(find(c2_right_led_idx), find(SIx));
c2_left_led_correct=intersect(find(c2_left_led_idx), find(SIx));
c2_right_led_correct_idx(c2_right_led_correct)=1;
c2_left_led_correct_idx(c2_left_led_correct)=1;
c3_right_led_correct=intersect(find(c3_right_led_idx), find(SIx));
c3_left_led_correct=intersect(find(c3_left_led_idx), find(SIx));
c3_right_led_correct_idx(c3_right_led_correct)=1;
c3_left_led_correct_idx(c3_left_led_correct)=1;
c4_right_led_correct=intersect(find(c4_right_led_idx), find(SIx));
c4_left_led_correct=intersect(find(c4_left_led_idx), find(SIx));
c4_right_led_correct_idx(c4_right_led_correct)=1;
c4_left_led_correct_idx(c4_left_led_correct)=1;
%% side sorted binofits
%control
[phat_c1_right_led, pci_c1_right_led] =binofit (sum(c1_right_led_correct_idx), length(c1_right_led));
[phat_c1_left_led, pci_c1_left_led] =binofit (sum(c1_left_led_correct_idx), length(c1_left_led));
[phat_c2_right_led, pci_c2_right_led] =binofit (sum(c2_right_led_correct_idx), length(c2_right_led));
[phat_c2_left_led, pci_c2_left_led] =binofit (sum(c2_left_led_correct_idx), length(c2_left_led));
[phat_c3_right_led, pci_c3_right_led] =binofit (sum(c3_right_led_correct_idx), length(c3_right_led));
[phat_c3_left_led, pci_c3_left_led] =binofit (sum(c3_left_led_correct_idx), length(c3_left_led));
[phat_c4_right_led, pci_c4_right_led] =binofit (sum(c4_right_led_correct_idx), length(c4_right_led));
[phat_c4_left_led, pci_c4_left_led] =binofit (sum(c4_left_led_correct_idx), length(c4_left_led));

[phat_c1_right_ctrl, pci_c1_right_ctrl]=binofit(sum(c1_right_ctrl_correct_idx), length(c1_right_ctrl));
[phat_c1_left_ctrl, pci_c1_left_ctrl]=binofit(sum(c1_left_ctrl_correct_idx), length(c1_left_ctrl));
[phat_c2_right_ctrl, pci_c2_right_ctrl]=binofit(sum(c2_right_ctrl_correct_idx), length(c2_right_ctrl));
[phat_c2_left_ctrl, pci_c2_left_ctrl]=binofit(sum(c2_left_ctrl_correct_idx), length(c2_left_ctrl));
[phat_c3_right_ctrl, pci_c3_right_ctrl]=binofit(sum(c3_right_ctrl_correct_idx), length(c3_right_ctrl));
[phat_c3_left_ctrl, pci_c3_left_ctrl]=binofit(sum(c3_left_ctrl_correct_idx), length(c3_left_ctrl));
[phat_c4_right_ctrl, pci_c4_right_ctrl]=binofit(sum(c4_right_ctrl_correct_idx), length(c4_right_ctrl));
[phat_c4_left_ctrl, pci_c4_left_ctrl]=binofit(sum(c4_left_ctrl_correct_idx), length(c4_left_ctrl));

right_ctrl_outcomes=[phat_c1_right_ctrl, phat_c2_right_ctrl, phat_c3_right_ctrl, phat_c4_right_ctrl];
left_ctrl_outcomes=[phat_c1_left_ctrl, phat_c2_left_ctrl, phat_c3_left_ctrl, phat_c4_left_ctrl];
right_led_outcomes=[phat_c1_right_led, phat_c2_right_led, phat_c3_right_led, phat_c4_right_led];
left_led_outcomes=[phat_c1_left_led, phat_c2_left_led, phat_c3_left_led, phat_c4_left_led];


%% calc pci's out
pci_c1_right_led=pci_c1_right_led(2)-pci_c1_right_led(1);
pci_c1_left_led=pci_c1_left_led(2)-pci_c1_left_led(1);
pci_c2_right_led=pci_c2_right_led(2)-pci_c2_right_led(1);
pci_c2_left_led=pci_c2_left_led(2)-pci_c2_left_led(1);
pci_c3_right_led=pci_c3_right_led(2)-pci_c3_right_led(1);
pci_c3_left_led=pci_c3_left_led(2)-pci_c3_left_led(1);
pci_c4_right_led=pci_c4_right_led(2)-pci_c4_right_led(1);
pci_c4_left_led=pci_c4_left_led(2)-pci_c4_left_led(1);

pci_c1_right_ctrl=pci_c1_right_ctrl(2)-pci_c1_right_ctrl(1);
pci_c1_left_ctrl=pci_c1_left_ctrl(2)-pci_c1_left_ctrl(1);
pci_c2_right_ctrl=pci_c2_right_ctrl(2)-pci_c2_right_ctrl(1);
pci_c2_left_ctrl=pci_c2_left_ctrl(2)-pci_c2_left_ctrl(1);
pci_c3_right_ctrl=pci_c3_right_ctrl(2)-pci_c3_right_ctrl(1);
pci_c3_left_ctrl=pci_c3_left_ctrl(2)-pci_c3_left_ctrl(1);
pci_c4_right_ctrl=pci_c4_right_ctrl(2)-pci_c4_right_ctrl(1);
pci_c4_left_ctrl=pci_c4_left_ctrl(2)-pci_c4_left_ctrl(1);


%%
% ******* beautiful bar graph here: 
bargroup1=[phat_c1_right_ctrl, phat_c1_left_ctrl, phat_c1_right_led, phat_c1_left_led];
bargroup2=[phat_c2_right_ctrl, phat_c2_left_ctrl, phat_c2_right_led, phat_c2_left_led];
bargroup3=[phat_c3_right_ctrl, phat_c3_left_ctrl, phat_c3_right_led, phat_c3_left_led];
bargroup4=[phat_c4_right_ctrl, phat_c4_left_ctrl, phat_c4_right_led, phat_c4_left_led];
fullgroups=[bargroup1; bargroup2; bargroup3; bargroup4];
errorgroup1=[phat_c1_right_ctrl, phat_c2_right_ctrl, phat_c3_right_ctrl, phat_c4_right_ctrl];
errorgroup2=[phat_c1_left_ctrl, phat_c2_left_ctrl, phat_c3_left_ctrl, phat_c4_left_ctrl];
errorgroup3=[phat_c1_right_led, phat_c2_right_led, phat_c3_right_led, phat_c4_right_led];
errorgroup4=[phat_c1_left_led, phat_c2_left_led, phat_c3_left_led, phat_c4_left_led];
x=1:4;
barwidth=1;
figure
hold on
b=bar(x,fullgroups,barwidth);
legend('Control Right', 'Control Left', 'LED Right', 'LED Left')
xticks([1 2 3 4]);
xticklabels({'101', '11', '2', '1.1'})
xtickangle(45)
xlabel('Contrast Ratio (Target/Distractor)')
ylabel('Percent chance of successful trial')
titlevar=[length(tGratingContrast), length(s_cropped_include)];
titlestr= sprintf(' Percent success rate of contrast discrimination task by side and contrast ratio in control and LED conditions. \nGenotype: SOM \nN = %d trials across %d sessions',titlevar(1), titlevar(2) ) 
title(titlestr);
b(1).FaceColor = [ 0.1 0.7 0.1];
b(2).FaceColor=[0 0.5 1];
b(3).FaceColor='red';
b(4).FaceColor=[1 1 0];
errors_ctrl_right=[pci_c1_right_ctrl, pci_c2_right_ctrl, pci_c3_right_ctrl, pci_c4_right_ctrl];
errors_ctrl_left=[pci_c1_left_ctrl, pci_c2_left_ctrl, pci_c3_left_ctrl, pci_c4_left_ctrl];
errors_led_right=[pci_c1_right_led, pci_c2_right_led, pci_c3_right_led, pci_c4_right_led];
errors_led_left=[pci_c1_left_led, pci_c2_left_led, pci_c3_left_led, pci_c4_left_led];
errorbar(x-0.275*barwidth,errorgroup1,errors_ctrl_right,'.k');
errorbar(x-0.09*barwidth,errorgroup2, errors_ctrl_left, '.k');
errorbar(x+0.09*barwidth,errorgroup3, errors_led_right, '.k');
errorbar(x+0.275*barwidth,errorgroup4,errors_led_left,'.k');





%% use binofit maybe get the same thing?

%get correct controls sum and length(correctcontrols), same for LED
%condition

correct_controls=intersect(find(control_idx), find(SIx)); %control iso
correct_controls_idx=zeros(1,length(tGratingContrast));
correct_controls_idx(correct_controls)=1;
correct_leds=intersect(find(tBlock2), find(SIx)); %LED iso
correct_leds_idx=zeros(1,length(tGratingContrast));
correct_leds_idx(correct_leds)=1;

correct_controls_c1_idx=zeros(1,length(tGratingContrast));
correct_controls_c2_idx=correct_controls_c1_idx;
correct_controls_c3_idx=correct_controls_c1_idx;
correct_controls_c4_idx=correct_controls_c1_idx;

correct_led_c1_idx=zeros(1,length(tGratingContrast));
correct_led_c2_idx=correct_led_c1_idx;
correct_led_c3_idx=correct_led_c1_idx;
correct_led_c4_idx=correct_led_c1_idx;

correct_controls_c1=intersect(find(correct_controls_idx), find(contrastratio==unqcontrast(4)));
correct_controls_c1_idx(correct_controls_c1)=1;
correct_controls_c2=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(3)));
correct_controls_c2_idx(correct_controls_c2)=1;
correct_controls_c3=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(2)));
correct_controls_c3_idx(correct_controls_c3)=1;
correct_controls_c4=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(1)));
correct_controls_c4_idx(correct_controls_c4)=1;



correct_led_c1=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(4)));
correct_led_c1_idx(correct_led_c1)=1;
correct_led_c2=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(3)));
correct_led_c2_idx(correct_led_c2)=1;
correct_led_c3=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(2)));
correct_led_c3_idx(correct_led_c3)=1;
correct_led_c4=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(1)));
correct_led_c4_idx(correct_led_c4)=1;

numcontrol_c1_trials=intersect(find(contrastratio==unqcontrast(4)), find(control_idx));
numcontrol_c2_trials=intersect(find(contrastratio==unqcontrast(3)), find(control_idx));
numcontrol_c3_trials=intersect(find(contrastratio==unqcontrast(2)), find(control_idx));
numcontrol_c4_trials=intersect(find(contrastratio==unqcontrast(1)), find(control_idx));

numled_c1_trials=intersect(find(tBlock2), find(contrastratio==unqcontrast(4)));
numled_c2_trials=intersect(find(tBlock2), find(contrastratio==unqcontrast(3)));
numled_c3_trials=intersect(find(tBlock2), find(contrastratio==unqcontrast(2)));
numled_c4_trials=intersect(find(tBlock2), find(contrastratio==unqcontrast(1)));

%%

[phat_controlc1,pci_controlc1]=binofit(sum(correct_controls_c1_idx), length(numcontrol_c1_trials));
[phat_controlc2, pci_controlc2]=binofit(sum(correct_controls_c2_idx), length(numcontrol_c2_trials));
[phat_controlc3,pci_controlc3]=binofit(sum(correct_controls_c3_idx), length(numcontrol_c3_trials));
[phat_controlc4,pci_controlc4]=binofit(sum(correct_controls_c4_idx), length(numcontrol_c4_trials));

%LED outcomes 
[phat_ledc1,pci_ledc1]=binofit(sum(correct_led_c1_idx), length(numled_c1_trials));
[phat_ledc2,pci_ledc2]=binofit(sum(correct_led_c2_idx), length(numled_c2_trials));
[phat_ledc3,pci_ledc3]=binofit(sum(correct_led_c3_idx), length(numled_c3_trials));
[phat_ledc4,pci_ledc4]=binofit(sum(correct_led_c4_idx), length(numled_c4_trials));


%led_outcome=[phat_ledc1, phat_ledc2, phat_ledc3, phat_ledc4];
%barci_control=[pci_controlc1(2)-pci_controlc1(1), pci_controlc2(2)-pci_controlc2(1), pci_controlc3(2)-pci_controlc3(1), pci_controlc4(2)-pci_controlc4(1)];
%barci_led=[pci_ledc1(2)-pci_ledc1(1), pci_ledc2(2)-pci_ledc2(1), pci_ledc3(2)-pci_ledc3(1), pci_ledc4(2)-pci_ledc4(1)];
%barci=[barci_control(1), barci_led(1), barci_control(2), barci_led(2), barci_control(3), barci_led(3), barci_control(4), barci_led(4)];
%%
data_to_plot=[control_outcome(1), led_outcome(1); 
    control_outcome(2), led_outcome(2); 
    control_outcome(3), led_outcome(3); 
    control_outcome(4), led_outcome(4)];


figure
hold on
x=1:4;
barwidth=1;
b=bar(x,data_to_plot,barwidth);
b(1).FaceColor = 'yellow';
b(2).FaceColor='cyan';
xtickangle(45)
xlabel('Contrast Ratio (Target/Distractor)')
ylabel('Percent chance of successful trial')
legend('Control', 'LED')
errorbar(x-0.15*barwidth,control_outcome,barci_control,'.k');
errorbar(x+0.15*barwidth,led_outcome,barci_led,'.k');


titlevar=[length(tGratingContrast), length(s_cropped_include)];
titlestr=sprintf('Percent success in Control and LED Conditions by Contrast Ratio. \nSubject: i414 \nGenotype: SOM \nLED Power= mW \nN = %d trials across %d sessions',titlevar(1), titlevar(2));
title(titlestr);
%set(gca,'xticklabel',labels)
%% APPENDIX - Unneeded but possibly useful code
%Appendix A
%{
 %create contrast sorted LED indices
c1_led_right_idx=zeros(1,length(tGratingContrast));
c1_led_left_idx=c1_led_right_idx;
c2_led_right_idx=c1_led_right_idx;
c2_led_left_idx=c1_led_right_idx;
c3_led_right_idx=c1_led_right_idx;
c3_led_left_idx=c1_led_right_idx;
c4_led_right_idx=c1_led_right_idx;
c4_led_left_idx=c1_led_right_idx;

c1_led_right=intersect(find(led_idx), find(contrastratio==c1));
%c1_led_left=intersect(find(led_idx

%}

% Appendix B
%{
for i=1:100 %bootstrap control samples only, randomly index 500 of them, then proceed
    random_idx=datasample(controls,500);
    sample_idx=sort(random_idx);
    success_sample_idx=SIx(sample_idx);
    sum_sample_success=sum(success_sample_idx);
    percent_sample(i)=sum_sample_success/length(random_idx);
end

for i=1:100
    random_idx2=datasample(leds,500);
 
    led_sample_idx=sort(random_idx2);
    led_success_idx=SIx(led_sample_idx);
    sum_led_success=sum(led_success_idx);
    percent_led_sample(i)=sum_led_success/length(random_idx2);
end

%} 
% Appendix C
%{
set1=[percent_sample]'; 
set2=[percent_led_sample]';
group=[ones(size(set1)); 2 * ones(size(set2))];

figure(1)
hold on
boxplot([set1; set2],group,'Notch', 'on')
ylim([0.65,0.95])
%}

% Appendix D

%histograms
%{
figure(2)
hold on
histogram(percent_sample)
histogram(percent_led_sample)
legend('Control', 'Optogenetic Stimulation')
title('Distribution of percent correct by trial condition in 100x bootstrap')
xlabel('Percent trials correct')
ylabel('Count')
are_means_diff=ttest2(percent_sample, percent_led_sample,'Vartype','unequal')
%}