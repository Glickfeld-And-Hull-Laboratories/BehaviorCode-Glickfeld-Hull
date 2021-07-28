%% all inclusive opto contrast analysis 

targetContrast=celleqel2mat_padded(input.tGratingContrast);
distractContrast = celleqel2mat_padded(input.dGratingContrast);
islefttrial=celleqel2mat_padded(input.tLeftTrial);
isrighttrial=~islefttrial;
leftresponse=celleqel2mat_padded(input.tLeftResponse);
rightresponse=celleqel2mat_padded(input.tRightResponse); 
ledtrials=celleqel2mat_padded(input.tBlock2TrialNumber);

contrastratio=targetContrast./distractContrast;
contrastratio=round(contrastratio,3);
unqcontrasts = unique(contrastratio);
c1 = unqcontrasts(1);
c2= unqcontrasts(2);
c3=unqcontrasts(3);
c4=unqcontrasts(4);

%% isolate non led trials and get #/% correct
regtrials=zeros(1, length(targetContrast));
regtrials(~ledtrials)=1;
%righttrials=zeros(1,length(targetContrast)); 
reglefttrialsidx=zeros(1,length(targetContrast));
regrighttrialsidx=zeros(1,length(targetContrast));
regrighttrials=intersect(find(regtrials), find(isrighttrial));
reglefttrials=intersect(find(regtrials), find(islefttrial));
reglefttrialsidx(reglefttrials)=1;
regrighttrialsidx(regrighttrials)=1;
ledlefttrials=intersect(find(ledtrials), find(islefttrial));
ledrighttrials=intersect(find(ledtrials), find(isrighttrial));
ledlefttrialsidx=zeros(1,length(targetContrast));
ledrighttrialsidx=ledlefttrialsidx;
ledlefttrialsidx(ledlefttrials)=1;
ledrighttrialsidx(ledrighttrials)=1;
ledrighttrialscorrect=intersect(find(ledrighttrialsidx), find(rightresponse));
ledlefttrialscorrect=intersect(find(ledlefttrialsidx), find(leftresponse));
ledleftcorrectidx=zeros(1,length(targetContrast));
ledrightcorrectidx=zeros(1,length(targetContrast));
ledleftcorrectidx(ledlefttrialscorrect)=1;
ledrightcorrectidx(ledrighttrialscorrect)=1;
%need to sort all led trials and correct ones by contrast

%%
ledrightc1=length(intersect(find(ledrighttrialsidx), find(contrastratio==c1)));
ledleftc1=length(intersect(find(ledlefttrialsidx), find(contrastratio==c1)));
ledrightc2=length(intersect(find(ledrighttrialsidx), find(contrastratio==c2)));
ledleftc2=length(intersect(find(ledlefttrialsidx), find(contrastratio==c2)));
ledrightc3=length(intersect(find(ledrighttrialsidx), find(contrastratio==c3)));
ledleftc3=length(intersect(find(ledlefttrialsidx), find(contrastratio==c3)));
ledrightc4=length(intersect(find(ledrighttrialsidx), find(contrastratio==c4)));
ledleftc4=length(intersect(find(ledlefttrialsidx), find(contrastratio==c4)));

ledcorrectrightc1=length(intersect(find(ledrightcorrectidx), find(contrastratio==c1)));
ledcorrectleftc1=length(intersect(find(ledleftcorrectidx), find(contrastratio==c1)));
ledcorrectrightc2=length(intersect(find(ledrightcorrectidx), find(contrastratio==c2)));
ledcorrectleftc2=length(intersect(find(ledleftcorrectidx), find(contrastratio==c2)));
ledcorrectrightc3=length(intersect(find(ledrightcorrectidx), find(contrastratio==c3)));
ledcorrectleftc3=length(intersect(find(ledleftcorrectidx), find(contrastratio==c3)));
ledcorrectrightc4=length(intersect(find(ledrightcorrectidx), find(contrastratio==c4)));
ledcorrectleftc4=length(intersect(find(ledleftcorrectidx), find(contrastratio==c4)));

ledpctrightc1=ledcorrectrightc1/ledrightc1; % led percent corrects
ledpctleftc1=ledcorrectleftc1/ledleftc1;
ledpctrightc2=ledcorrectrightc2/ledrightc2;
ledpctleftc2=ledcorrectleftc2/ledleftc2;
ledpctrightc3=ledcorrectrightc3/ledrightc3;
ledpctleftc3=ledcorrectleftc3/ledleftc3;
ledpctrightc4=ledcorrectrightc4/ledrightc4;
ledpctleftc4=ledcorrectleftc4/ledleftc4;

ledpctsarray=[ledpctleftc4, ledpctleftc3, ledpctleftc2, ledpctleftc1, ledpctrightc1, ledpctrightc2, ledpctrightc3, ledpctrightc4];



%% get correct reg trials and % 

%% sort percents and contrasts into arrays

%% plot
plotcontrasts=[ 1/c4, 1/c3, 1/c2, 1/c1, c1, c2,c3,c4];
plotcontrasts=round(plotcontrasts, 3);
numbins=length(plotcontrasts);
pc_prime=[1:numbins];
%pctvalues=smooth(pctvalues); % unsure if needed or valuable? 
figure(1)
plot(pc_prime , ledpctsarray)
%plot(pc_prime, pctvalues)%'-o', 'MarkerFaceColor', [0,0,0.5])
set(gca,'XtickLabel', plotcontrasts);

