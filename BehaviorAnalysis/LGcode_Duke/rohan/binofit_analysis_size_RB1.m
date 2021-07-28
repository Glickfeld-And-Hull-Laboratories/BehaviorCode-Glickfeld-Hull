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




%% CONTRAST SORTING

sizeratio=tGratingDiameter./dGratingDiameter;
contrastratio=round(sizeratio,4);
contrast1_rightidx=zeros(1, length(tGratingDiameter)); %declare contrast idx arrays
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
%% use binofit maybe get the same thing?

%get correct controls sum and length(correctcontrols), same for LED
%condition

correct_controls=intersect(find(control_idx), find(SIx)); %control iso
correct_controls_idx=zeros(1,length(tGratingDiameter));
correct_controls_idx(correct_controls)=1;
correct_leds=intersect(find(tBlock2), find(SIx)); %LED iso
correct_leds_idx=zeros(1,length(tGratingDiameter));
correct_leds_idx(correct_leds)=1;

correct_controls_c1_idx=zeros(1,length(tGratingDiameter));
correct_controls_c2_idx=correct_controls_c1_idx;
correct_controls_c3_idx=correct_controls_c1_idx;
correct_controls_c4_idx=correct_controls_c1_idx;

correct_led_c1_idx=zeros(1,length(tGratingDiameter));
correct_led_c2_idx=correct_led_c1_idx;
correct_led_c3_idx=correct_led_c1_idx;
correct_led_c4_idx=correct_led_c1_idx;
correct_ctrl_c1_right_idx=correct_led_c1_idx;
correct_ctrl_c1_left_idx=correct_led_c1_idx;
correct_ctrl_c2_right_idx=correct_led_c1_idx;
correct_ctrl_c2_left_idx=correct_led_c1_idx;
correct_ctrl_c3_right_idx=correct_led_c1_idx;
correct_ctrl_c3_left_idx=correct_led_c1_idx;
correct_ctrl_c4_right_idx=correct_led_c1_idx;
correct_ctrl_c4_left_idx=correct_led_c1_idx;



correct_controls_c1=intersect(find(correct_controls_idx), find(contrastratio==unqcontrast(4)));
correct_controls_c1_idx(correct_controls_c1)=1;
correct_controls_c2=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(3)));
correct_controls_c2_idx(correct_controls_c2)=1;
correct_controls_c3=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(2)));
correct_controls_c3_idx(correct_controls_c3)=1;
correct_controls_c4=intersect(find(correct_controls_idx),find(contrastratio==unqcontrast(1)));
correct_controls_c4_idx(correct_controls_c4)=1;

correct_ctrls_c1_right=intersect(find(correct_controls_c1_idx), find(~tLeftTrial));
correct_ctrls_c1_left=intersect(find(correct_controls_c1_idx), find(tLeftTrial));
correct_ctrls_c1_right_idx(correct_ctrls_c1_right)=1;
correct_ctrls_c1_left_idx(correct_ctrls_c1_left)=1;
correct_ctrls_c2_right=intersect(find(correct_controls_c2_idx), find(~tLeftTrial));
correct_ctrls_c2_left=intersect(find(correct_controls_c2_idx), find(tLeftTrial));
correct_ctrls_c2_right_idx(correct_ctrls_c2_right)=1;
correct_ctrls_c2_left_idx(correct_ctrls_c2_left)=1;
correct_ctrls_c3_right=intersect(find(correct_controls_c3_idx), find(~tLeftTrial));
correct_ctrls_c3_left=intersect(find(correct_controls_c3_idx), find(tLeftTrial));
correct_ctrls_c3_right_idx(correct_ctrls_c3_right)=1;
correct_ctrls_c3_left_idx(correct_ctrls_c3_left)=1;
correct_ctrls_c4_right=intersect(find(correct_controls_c4_idx), find(~tLeftTrial));
correct_ctrls_c4_left=intersect(find(correct_controls_c4_idx), find(tLeftTrial));
correct_ctrls_c4_right_idx(correct_ctrls_c4_right)=1;
correct_ctrls_c4_left_idx(correct_ctrls_c4_left)=1;




%%

correct_led_c1=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(4)));
correct_led_c1_idx(correct_led_c1)=1;
correct_led_c2=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(3)));
correct_led_c2_idx(correct_led_c2)=1;
correct_led_c3=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(2)));
correct_led_c3_idx(correct_led_c3)=1;
correct_led_c4=intersect(find(correct_leds_idx), find(contrastratio==unqcontrast(1)));
correct_led_c4_idx(correct_led_c4)=1;

%% maybe need? 
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
control_outcome_c1= [phat_controlc1 pci_controlc1];
[phat_controlc2, pci_controlc2]=binofit(sum(correct_controls_c2_idx), length(numcontrol_c2_trials));
control_outcome_c2=[phat_controlc2, pci_controlc2];
[phat_controlc3,pci_controlc3]=binofit(sum(correct_controls_c3_idx), length(numcontrol_c3_trials));
control_outcome_c3=[phat_controlc3,pci_controlc3];
[phat_controlc4,pci_controlc4]=binofit(sum(correct_controls_c4_idx), length(numcontrol_c4_trials));
control_outcome_c4=[phat_controlc4,pci_controlc4];
%LED outcomes below
[phat_ledc1,pci_ledc1]=binofit(sum(correct_led_c1_idx), length(numled_c1_trials));
led_outcome_c1=[phat_ledc1,pci_ledc1];
[phat_ledc2,pci_ledc2]=binofit(sum(correct_led_c2_idx), length(numled_c2_trials));
led_outcome_c2=[phat_ledc2,pci_ledc2];
[phat_ledc3,pci_ledc3]=binofit(sum(correct_led_c3_idx), length(numled_c3_trials));
led_outcome_c3=[phat_ledc3,pci_ledc3];
[phat_ledc4,pci_ledc4]=binofit(sum(correct_led_c4_idx), length(numled_c4_trials));
led_outcome_c4=[phat_ledc4,pci_ledc4];

control_outcome=[phat_controlc1, phat_controlc2, phat_controlc3, phat_controlc4];
led_outcome=[phat_ledc1, phat_ledc2, phat_ledc3, phat_ledc4];
barci_control=[pci_controlc1(2)-pci_controlc1(1), pci_controlc2(2)-pci_controlc2(1), pci_controlc3(2)-pci_controlc3(1), pci_controlc4(2)-pci_controlc4(1)];
barci_led=[pci_ledc1(2)-pci_ledc1(1), pci_ledc2(2)-pci_ledc2(1), pci_ledc3(2)-pci_ledc3(1), pci_ledc4(2)-pci_ledc4(1)];
barci=[barci_control(1), barci_led(1), barci_control(2), barci_led(2), barci_control(3), barci_led(3), barci_control(4), barci_led(4)];
%%
data_to_plot=[control_outcome(1), led_outcome(1); 
    control_outcome(2), led_outcome(2); 
    control_outcome(3), led_outcome(3); 
    control_outcome(4), led_outcome(4)];


figure
hold on
x=1:4;
barwidth=1
b=bar(x,data_to_plot,barwidth);
b(1).FaceColor = 'yellow';
b(2).FaceColor='cyan';
xtickangle(45)
xlabel('Size Ratio (Target/Distractor)')
ylabel('Percent chance of successful trial')
legend('Control', 'LED')
errorbar(x-0.15*barwidth,control_outcome,barci_control,'.k');
errorbar(x+0.15*barwidth,led_outcome,barci_led,'.k');


titlevar=[length(tGratingContrast), length(s_cropped_include)];
titlestr=sprintf('Binomial Estimate of Success Probability in Control and LED Conditions by Size Ratio. \nSubject: i600 \nGenotype: PV \nN = %d trials across %d sessions',titlevar(1), titlevar(2));
title(titlestr);
%set(gca,'xticklabel',labels)


