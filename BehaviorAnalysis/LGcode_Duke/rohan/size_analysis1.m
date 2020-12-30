%%
uniquesizes=unique(sizes)

%% copied from contrast analysis for size analysis use
islefttrial=celleqel2mat_padded(input.tLeftTrial)
distractSize=celleqel2mat_padded(input.dSize)

nums1rightcorrect=length(intersect(find(s1righttrials), find(rightresponse)));
pcts1correctright=nums1rightcorrect/length(s1righttrials);
nums1leftcorrect=length(intersect(find(s1lefttrials), find(leftresponse))); 
pcts1correctleft=nums1leftcorrect/length(s1lefttrials);
nums2rightcorrect=length(intersect(find(s2righttrials), find(rightresponse)));
pcts2correctright=nums2rightcorrect/length(s2righttrials);
nums2leftcorrect=length(intersect(find(s2lefttrials), find(leftresponse)));
pcts2correctleft=nums2leftcorrect/length(s2lefttrials);
nums3rightcorrect=length(intersect(find(s3righttrials), find(rightresponse)));
pcts3correctright=nums3rightcorrect/length(s3righttrials);
nums3leftcorrect=length(intersect(find(s3lefttrials), find(leftresponse)));
pcts3correctleft=nums3leftcorrect/length(s3lefttrials);
nums4rightcorrect=length(intersect(find(s4righttrials), find(rightresponse)));
pctc4correctright=numc4rightcorrect/length(s4righttrials);
nums4leftcorrect=length(intersect(find(s4lefttrials), find(leftresponse)));
pcts4correctleft=nums4leftcorrect/length(s4lefttrials);

pctvalues=[ pcts4correctleft, pcts3correctleft,pcts2correctleft, pcts1correctleft, pcts1correctright, pcts2correctright, pcts3correctright, pcts4correctright];

plotsizes=[ 1/s4, 1/s3, 1/s2, 1/s1, s1, s2,s3,s4];
plotsizess=round(plotsizes, 3);
numbins=length(plotsizes);
ps_prime=[1:numbins];
%pctvalues=smooth(pctvalues); % unsure if needed or valuable? 
figure(1)
plot(ps_prime, pctvalues)%'-o', 'MarkerFaceColor', [0,0,0.5])
set(gca,'XtickLabel', plotsizes);
%xlim([0, 105])


