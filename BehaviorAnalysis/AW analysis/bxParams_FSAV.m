%% cutoffs and binning
early_cutoff = 0.35;
lapse_cutoff = 0.9;
minTrN_ms = 15;
minTrN_expt = 2;
nBins = 6;
minCyclesFA = 2;
minTrialLengthMs = 500;

visThreshAll = 30;
highThreshold = 0.75;
audThreshAll = 0.015;

timeBins = [0 1100 5000];

%% structure ID
visualTrials = 1;
auditoryTrials = 2;
valid  = 1;
invalid = 2;

%% colors
cueColor = {[0 0 0];[.5 .5 1]};
AVColor = {[0 0 0];[0.5 0.5 0.5]};
hiLoColor = {[0.5 0.5 0.5];[0 0 0]};

%% significance testing
attnTestAlpha = 0.05;