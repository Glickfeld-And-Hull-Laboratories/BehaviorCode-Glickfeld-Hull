%attempt to plot percent correct at each target ratio by each side
% does weber's law hold at lower contrast ratios? maybe not - greg

%t=target d=distractor
%% load data into arrays
%celleqel2mat_padded(input. );

%% try to sort percent correct by contrast ratio ****

targetContrast=celleqel2mat_padded(input.tGratingContrast);
distractContrast = celleqel2mat_padded(input.dGratingContrast);
islefttrial=celleqel2mat_padded(input.tLeftTrial);
isrighttrial=~islefttrial;
leftresponse=celleqel2mat_padded(input.tLeftResponse);
rightresponse=celleqel2mat_padded(input.tRightResponse);

contrastratio=targetContrast./distractContrast;

unqcontrasts = unique(contrastratio);
%unqcontrasts= round(unqcontrasts,3);
c1 = unqcontrasts(1);
c2= unqcontrasts(2);
c3=unqcontrasts(3);
c4=unqcontrasts(4);

corrects=zeros(length(islefttrial),1);

correctleft=intersect(find(islefttrial), find(leftresponse));
correctright=intersect(find(isrighttrial), find(rightresponse)); 
numcorrectleft=length(correctleft);
numcorrectright=length(correctright);
corrects(correctleft)=1;
corrects(correctright)=1;

c1correct=intersect(find(contrastratio==c1), find(corrects));
numc1correct=length(c1correct);
totalc1trials=length(intersect(find(contrastratio==c1), find(targetContrast)));
percentc1correct=numc1correct/totalc1trials;

c2correct=intersect(find(contrastratio==c2), find(corrects));
numc2correct=length(c2correct);
totalc2trials=length(intersect(find(contrastratio==c2), find(targetContrast)));
percentc2correct=numc2correct/totalc2trials;

c3correct=intersect(find(contrastratio==c3), find(corrects));
numc3correct=length(c3correct);
totalc3trials=length(intersect(find(contrastratio==c3), find(targetContrast)));%percentc3correct=numc3correct/totalc3trials;
percentc3correct=numc3correct/totalc3trials;

c4correct=intersect(find(contrastratio==c4), find(corrects));
numc4correct=length(c4correct);
totalc4trials=length(intersect(find(contrastratio==c4), find(targetContrast)));
percentc4correct=numc4correct/totalc4trials;
%%
figure; 
plot(unqcontrasts, [percentc1correct percentc2correct percentc3correct percentc4correct],'-o') 


%% plot percent correct by right left and contrasts (similar to curve in PDF) 

c1righttrials=intersect(find(contrastratio==c1), find(isrighttrial));
c1lefttrials=intersect(find(contrastratio==c1), find(islefttrial)); 
c2righttrials=intersect(find(contrastratio==c2), find(isrighttrial)); 
c2lefttrials=intersect(find(contrastratio==c2), find(islefttrial));
c3righttrials=intersect(find(contrastratio==c3), find(isrighttrial));
c3lefttrials=intersect(find(contrastratio==c3), find(islefttrial));
c4righttrials=intersect(find(contrastratio==c4), find(isrighttrial));
c4lefttrials=intersect(find(contrastratio==c4), find(islefttrial));

actualc1right=zeros(length(targetContrast), 1); 
actualc1left=actualc1right;
actualc2right=actualc1right;
actualc2left=actualc1right;
actualc3right=actualc1right;
actualc3left=actualc1right;
actualc4right=actualc1right;
actualc4left=actualc1right;


actualc1right(c1righttrials)=1;
actualc1left(c1lefttrials)=1;
actualc2right(c2righttrials)=1;
actualc2left(c2lefttrials)=1;
actualc3right(c3righttrials)=1;
actualc3left(c3lefttrials)=1;
actualc4right(c4righttrials)=1;
actualc4left(c4lefttrials)=1;

ignoresidx=zeros(length(targetContrast), 1);
ignores= intersect(find(~rightresponse), find(~leftresponse));
ignoresidx(ignores)=1;
ignoresright=intersect(find(ignoresidx),find(isrighttrial));
ignoresleft=intersect(find(ignoresidx), find(islefttrial));
ignoresrightidx=zeros(length(targetContrast), 1);
ignoresleftidx=zeros(length(targetContrast), 1);
ignoresrightidx(ignoresright)=1;
ignoresleftidx(ignoresleft)=1;


numc1rightcorrect=length(intersect(find(c1righttrials), find(rightresponse)));
c1rightignores = intersect(find(contrastratio==c1), find(ignoresrightidx));
pctc1correctright=numc1rightcorrect/(length(c1righttrials)-length(c1rightignores));
numc1leftcorrect=length(intersect(find(c1lefttrials), find(leftresponse))); 
c1leftignores=intersect(find(contrastratio==c1), find(ignoresleftidx));
pctc1correctleft=numc1leftcorrect/(length(c1lefttrials)-length(c1leftignores));
numc2rightcorrect=length(intersect(find(c2righttrials), find(rightresponse)));
c2rightignores=intersect(find(contrastratio==c2), find(ignoresrightidx));
pctc2correctright=numc2rightcorrect/(length(c2righttrials)-length(c2rightignores));
numc2leftcorrect=length(intersect(find(c2lefttrials), find(leftresponse)));
c2leftignores=intersect(find(contrastratio==c2), find(ignoresleftidx));
pctc2correctleft=numc2leftcorrect/(length(c2lefttrials)-length(c2leftignores));
numc3rightcorrect=length(intersect(find(c3righttrials), find(rightresponse)));
c3rightignores=intersect(find(contrastratio==c3), find(ignoresrightidx)); 
pctc3correctright=numc3rightcorrect/(length(c3righttrials)-length(c3rightignores));
numc3leftcorrect=length(intersect(find(c3lefttrials), find(leftresponse)));
c3leftignores=intersect(find(contrastratio==c3), find(ignoresleftidx));
pctc3correctleft=numc3leftcorrect/(length(c3lefttrials)-length(c3leftignores));
numc4rightcorrect=length(intersect(find(c4righttrials), find(rightresponse)));
c4rightignores=intersect(find(contrastratio==c4), find(ignoresrightidx));
pctc4correctright=numc4rightcorrect/(length(c4righttrials)-length(c4rightignores));
numc4leftcorrect=length(intersect(find(c4lefttrials), find(leftresponse)));
c4leftignores=intersect(find(contrastratio==c4), find(ignoresleftidx));
pctc4correctleft=numc4leftcorrect/(length(c4lefttrials)-length(c4leftignores));

%% plotting
pctvalues=[ pctc4correctleft, pctc3correctleft,pctc2correctleft, pctc1correctleft, pctc1correctright, pctc2correctright, pctc3correctright, pctc4correctright];
plotcontrasts=[ 1/c4, 1/c3, 1/c2, 1/c1, c1, c2,c3,c4];
plotcontrasts=round(plotcontrasts, 3);
numbins=length(plotcontrasts);
pc_prime=[1:numbins];
%pctvalues=smooth(pctvalues); % unsure if needed or valuable? 
figure(1)
plot(pc_prime, pctvalues)%'-o', 'MarkerFaceColor', [0,0,0.5])
set(gca,'XtickLabel', plotcontrasts);
%xlim([0, 105])











