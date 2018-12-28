%% cutoffs and binning
early_cutoff = 0.35;
lapse_cutoff = 0.9;
minTrN_ms = 20;
minTrN_expt = 2;
nBins = 6;
minCyclesFA = 2;
minTrialLengthMs = 500;
visRTwindow = [200 550];
audRTwindow = [150 450];

visThreshAll = 30;
highThreshold = 0.8;
audThreshAll = 0.015;

timeBins = [0 1100 5000];

minVisTarget = 8;
maxVisTarget = 90;
minAudTarget = 0.0003;
maxAudTarget = 1;
visBinEdges = exp(linspace(log(minVisTarget-1),log(maxVisTarget),nBins+1));
audBinEdges = exp(linspace(log(minAudTarget-(0.5*minAudTarget)),...
    log(maxAudTarget),nBins+2));
audBinEdges = audBinEdges([1,3:end]);
%% structure ID
visualTrials = 1;
auditoryTrials = 2;
valid  = 1;
invalid = 2;

%% colors & lines
cueColor = {[0 0 0];[.5 .5 1]};
AVColor = {[0 0 0];[0.5 0.5 0.5]};
hiLoColor = {[0.5 0.5 0.5];[0 0 0]};
rewardedLine = {'-';':'};
rewardedFaceColor = {[0 0 0];'none'};
%% significance testing
attnTestAlpha = 0.05;
RTanovaAlpha = 0.05;
sampleSizePower = 0.95;