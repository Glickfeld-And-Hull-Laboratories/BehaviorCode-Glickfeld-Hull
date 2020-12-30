
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

contrastratio=tGratingContrast./dGratingContrast;
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
c4_right_led=intersect(find(contrast4_rightidx), find(led_idx));
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

% SOME LEFT AND RIGHT LED PHAT's are equal, LOOK INTO IT

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
%%
bargroup1=[phat_c1_right_ctrl-phat_c1_right_led, phat_c1_left_ctrl-phat_c1_left_led];
bargroup2=[phat_c2_right_ctrl-phat_c2_right_led, phat_c2_left_ctrl-phat_c2_left_led];
bargroup3=[phat_c3_right_ctrl-phat_c3_right_led, phat_c3_left_ctrl-phat_c3_left_led];
bargroup4=[phat_c4_right_ctrl-phat_c4_right_led, phat_c4_left_ctrl-phat_c4_left_led];
fullgroups=[bargroup1; bargroup2; bargroup3; bargroup4];
%%
errorgroup1=[phat_c1_right_ctrl, phat_c2_right_ctrl, phat_c3_right_ctrl, phat_c4_right_ctrl];
errorgroup2=[phat_c1_left_ctrl, phat_c2_left_ctrl, phat_c3_left_ctrl, phat_c4_left_ctrl];
errorgroup3=[phat_c1_right_led, phat_c2_right_led, phat_c3_right_led, phat_c4_right_led];
errorgroup4=[phat_c1_left_led, phat_c2_left_led, phat_c3_left_led, phat_c4_left_led];

%% binocdf use section
binotest_c1_left=binocdf(length(c1_left_led_correct), length(c1_left_led), phat_c1_left_ctrl);
binotest_c2_left=binocdf(length(c2_left_led_correct), length(c2_left_led), phat_c2_left_ctrl);
binotest_c3_left=binocdf(length(c3_left_led_correct), length(c3_left_led), phat_c3_left_ctrl);
binotest_c4_left=binocdf(length(c4_left_led_correct), length(c4_left_led), phat_c4_left_ctrl);

binotest_c1_right=binocdf(length(c1_right_led_correct), length(c1_right_led), phat_c1_right_ctrl);
binotest_c2_right=binocdf(length(c2_right_led_correct), length(c2_right_led), phat_c2_right_ctrl);
binotest_c3_right=binocdf(length(c3_right_led_correct), length(c3_right_led), phat_c3_right_ctrl);
binotest_c4_right=binocdf(length(c4_right_led_correct), length(c4_right_led), phat_c4_right_ctrl);

%%
x=1:4;
barwidth=1;
figure
hold on
b=bar(x,fullgroups,barwidth);
legend('\Delta Right (Control-LED)', '\Delta Left (Control-LED)')
xticks([1 2 3 4]);
xticklabels({'101', '11', '2', '1.1'})
xtickangle(45)
b(1).FaceColor='black'
b(2).FaceColor='blue'
xlabel('Contrast Ratio (Target/Distractor)')
ylabel('Percent success')
ctr1 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
ctr2 = bsxfun(@plus, b(2).XData, [b(2).XOffset]');
ctrs=horzcat(ctr1, ctr2);
ctrs=sort(ctrs);
clear ctr1; clear ctr2;
if binotest_c1_right<0.05  | binotest_c1_right>0.95
  %plot(ctrs(1:2), [2 2]*fullgroups(1,2)*3, '-k', 'LineWidth',2)
    plot(ctrs(1), fullgroups(1,1)*1.1, '*k')
end
if binotest_c1_left<0.05 | binotest_c1_left>0.95
  %plot(ctrs(1:2), [2 2]*fullgroups(1,2)*3, '-k', 'LineWidth',2)
    plot(ctrs(2), fullgroups(1,2)*1.1, '*k')
end
if binotest_c2_right<0.05  | binotest_c2_right>0.95
 %  plot(ctrs(3:4), [2 2]*fullgroups(1,2)*3, '-k', 'LineWidth',2)
    plot(ctrs(3), fullgroups(2,1)*1.1, '*k')
end
if binotest_c2_left<0.05  | binotest_c2_left>0.95
    plot(ctrs(4), fullgroups(2,2)*1.1, '*k')
end
if binotest_c3_right<0.05  | binotest_c3_right>0.95
    plot(ctrs(5), fullgroups(3,1)*1.1, '*k')
end
if binotest_c3_left<0.05  | binotest_c3_left>0.95
    plot(ctrs(6), fullgroups(3,2)*1.1, '*k')
end
if binotest_c4_right<0.05  | binotest_c4_right>0.95
    plot(ctrs(7), fullgroups(4,1)*1.1, '*k')
end

if binotest_c4_left<0.05  | binotest_c4_left>0.95
    plot(ctrs(8), fullgroups(4,2)*1.1, '*k')
end


titlevar=[length(tGratingContrast), length(s_cropped_include), s_cropped_include.laserpower];
titlestr= sprintf('i428 Difference across control and LED conditions in percent success rate \n of contrast discrimination by side and T/D contrast ratio. \nGenotype: PV \nLaser Power = %g mW \nN = %d trials across %d sessions',titlevar(3),titlevar(1), titlevar(2) ) 
title(titlestr);


%%
errors_ctrl_right=[pci_c1_right_ctrl, pci_c2_right_ctrl, pci_c3_right_ctrl, pci_c4_right_ctrl];
errors_ctrl_left=[pci_c1_left_ctrl, pci_c2_left_ctrl, pci_c3_left_ctrl, pci_c4_left_ctrl];
errors_led_right=[pci_c1_right_led, pci_c2_right_led, pci_c3_right_led, pci_c4_right_led];
errors_led_left=[pci_c1_left_led, pci_c2_left_led, pci_c3_left_led, pci_c4_left_led];
errorbar(x-0.275*barwidth,errorgroup1,errors_ctrl_right,'.k');
errorbar(x-0.09*barwidth,errorgroup2, errors_ctrl_left, '.k');
errorbar(x+0.09*barwidth,errorgroup3, errors_led_right, '.k');
errorbar(x+0.275*barwidth,errorgroup4,errors_led_left,'.k');


