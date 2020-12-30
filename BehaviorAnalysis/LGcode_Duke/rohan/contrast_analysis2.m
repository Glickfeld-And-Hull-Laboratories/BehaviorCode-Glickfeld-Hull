%% section 1 - inputs
targetContrast=celleqel2mat_padded(input.tGratingContrast);
distractContrast = celleqel2mat_padded(input.dGratingContrast);
islefttrial=celleqel2mat_padded(input.tLeftTrial);
isrighttrial=~islefttrial;
leftresponse=celleqel2mat_padded(input.tLeftResponse);
rightresponse=celleqel2mat_padded(input.tRightResponse); 
ledtrials=celleqel2mat_padded(input.tBlock2TrialNumber);

contrastratio=targetContrast./distractContrast;

unqcontrasts = unique(contrastratio);
%unqcontrasts= round(unqcontrasts,3);
c1 = unqcontrasts(1);
c2= unqcontrasts(2);
c3=unqcontrasts(3);
c4=unqcontrasts(4);


%%  Section 2 - analyze
%plot percent correct by right left and contrasts (similar to curve in PDF) 

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


numc1rightcorrect=length(intersect(c1righttrials, find(rightresponse)));
c1rightignores = intersect(find(contrastratio==c1), find(ignoresrightidx));
pctc1correctright=numc1rightcorrect/(length(c1righttrials)-length(c1rightignores));
numc1leftcorrect=length(intersect(c1lefttrials, find(leftresponse))); 
c1leftignores=intersect(find(contrastratio==c1), find(ignoresleftidx));
pctc1correctleft=numc1leftcorrect/(length(c1lefttrials)-length(c1leftignores));
numc2rightcorrect=length(intersect(c2righttrials, find(rightresponse)));
c2rightignores=intersect(find(contrastratio==c2), find(ignoresrightidx));
pctc2correctright=numc2rightcorrect/(length(c2righttrials)-length(c2rightignores));
numc2leftcorrect=length(intersect(c2lefttrials, find(leftresponse)));
c2leftignores=intersect(find(contrastratio==c2), find(ignoresleftidx));
pctc2correctleft=numc2leftcorrect/(length(c2lefttrials)-length(c2leftignores));
numc3rightcorrect=length(intersect(c3righttrials, find(rightresponse)));
c3rightignores=intersect(find(contrastratio==c3), find(ignoresrightidx)); 
pctc3correctright=numc3rightcorrect/(length(c3righttrials)-length(c3rightignores));
numc3leftcorrect=length(intersect(c3lefttrials, find(leftresponse)));
c3leftignores=intersect(find(contrastratio==c3), find(ignoresleftidx));
pctc3correctleft=numc3leftcorrect/(length(c3lefttrials)-length(c3leftignores));
numc4rightcorrect=length(intersect(c4righttrials, find(rightresponse)));
c4rightignores=intersect(find(contrastratio==c4), find(ignoresrightidx));
pctc4correctright=numc4rightcorrect/(length(c4righttrials)-length(c4rightignores));
numc4leftcorrect=length(intersect(c4lefttrials, find(leftresponse)));
c4leftignores=intersect(find(contrastratio==c4), find(ignoresleftidx));
pctc4correctleft=numc4leftcorrect/(length(c4lefttrials)-length(c4leftignores));

%% opto trials
c1ledrightidx=zeros(length(actualc1right), 1);
c2ledrightidx=zeros(length(actualc1right), 1);
c3ledrightidx=zeros(length(actualc1right), 1);
c4ledrightidx=zeros(length(actualc1right), 1);
c1ledleftidx=zeros(length(actualc1right), 1);
c2ledlefttidx=zeros(length(actualc1right), 1);
c3ledleftidx=zeros(length(actualc1right), 1);
c4ledleftidx=zeros(length(actualc1right), 1);

ledleftidx=zeros(length(actualc1right), 1);
ledrightidx=zeros(length(actualc1right),1);
ledlefttrials=intersect(find(ledtrials), find(islefttrial));
ledrighttrials=intersect(find(ledtrials), find(isrighttrial));
ledleftidx(ledlefttrials)=1;
ledrightidx(ledrighttrials)=1;
%%

c1ledright=intersect(find(contrastratio==c1), find(ledrightidx)); %%right led trials indexing by contrast
c1ledrightidx(c1ledright)=1;
c2ledright=intersect(find(contrastratio==c2), find(ledrightidx));
c2ledrightidx(c2ledright)=1;
c3ledright=intersect(find(contrastratio==c3), find(ledrightidx));
c3ledrightidx(c3ledright)=1;
c4ledright=intersect(find(contrastratio==c4), find(ledrightidx));
c4ledrightidx(c4ledright)=1;

c1ledleft=intersect(find(contrastratio==c1), find(ledleftidx)); %%left led trials indexing by contrast
c1ledleftidx(c1ledleft)=1;
c2ledleft=intersect(find(contrastratio==c2), find(ledleftidx));
c2ledleftidx(c2ledleft)=1;
c3ledleft=intersect(find(contrastratio==c3), find(ledleftidx));
c3ledleftidx(c3ledleft)=1;
c4ledleft=intersect(find(contrastratio==c4), find(ledleftidx));
c4ledleftidx(c4ledleft)=1;

numc1ledrightcorrect=length(intersect(find(c1ledrightidx), find(rightresponse)));
c1ledrightignores=intersect(find(c1ledright), find(ignoresrightidx));
pctc1ledright=numc1ledrightcorrect/(length(c1ledright)-length(c1ledrightignores));
numc2ledrightcorrect=length(intersect(find(c2ledrightidx), find(rightresponse)));
c2ledrightignores=intersect(find(c2ledright), find(ignoresrightidx));
%
%% plotting
pctvalues=[ pctc4correctleft, pctc3correctleft,pctc2correctleft, pctc1correctleft, pctc1correctright, pctc2correctright, pctc3correctright, pctc4correctright];
plotcontrasts=[ 1/c4, 1/c3, 1/c2, 1/c1, c1, c2,c3,c4];
plotcontrasts=round(plotcontrasts, 3);
numbins=length(plotcontrasts);
pc_prime=[1:numbins];
%pctvalues=smooth(pctvalues); % unsure if needed or valuable? 
figure(1)
plot(pc_prime, pctvalues,'-o', 'MarkerFaceColor', [0,0,0.5])
set(gca,'XtickLabel', plotcontrasts);
ylim([0,1]);
%xlim([0, 105])
xlabel('R/L Contrast'); ylabel('Percent trials correct');
title('Subject i414; performance on contrast discrimination task by contrast differences')