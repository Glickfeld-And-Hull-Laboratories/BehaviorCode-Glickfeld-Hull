clear all;close all
% awFSAVdatasets_longStimON_V1
% for slct_expt = 18:20
% ms = expt(slct_expt).SubNum;
% nrun = expt(slct_expt).nrun;
% exptInfo = cell(3,nrun);
% exptInfo(1,:) = repmat({expt(slct_expt).date},[1,nrun]);
% exptInfo(2,:) = cellstr(expt(slct_expt).time_mat)';

doSave = 0;
fnout = 'Z:\Analysis\_temp figs\180201';

ms = '756';
exptInfo = [{[];'1511';[]}];
nrun = size(exptInfo,2);
exptInfo(1,:) = repmat({'180821'},[1,nrun]);

for iexp = 1:size(exptInfo,2)
    dt = exptInfo{1,iexp};
    t = exptInfo{2,iexp};

d = loadMworksFile(ms,dt,t);

%%
trialRange = exptInfo{3,iexp}; % use [] if all trials

visTargets = celleqel2mat_padded(d.tGratingDirectionDeg);
audTargets = celleqel2mat_padded(d.tSoundTargetAmplitude);

invTrialInd = cell2mat(d.tShortCatchTrial);

invVisTargets = celleqel2mat_padded(d.tCatchGratingDirectionDeg);
invAudTargets = celleqel2mat_padded(d.tSoundCatchAmplitude);

valOutcome = d.trialOutcomeCell;
if ~isfield(d,'catchTrialOutcomeCell')
    invOutcome = nan(1,length(invVisTargets));
else
    invOutcome = d.catchTrialOutcomeCell;
end

hit = strcmp(valOutcome,'success');
miss = strcmp(valOutcome,'ignore');
fa = strcmp(valOutcome,'failure');

allVisTargets = unique(visTargets);
allAudTargets = unique(audTargets);

maxTargets = ismember(visTargets,allVisTargets(end-1:end)) | ismember(audTargets,allAudTargets(end-1:end));
maxTargetsTrialN = find(maxTargets);

figure
subplot 221
plot(maxTargetsTrialN,hit(maxTargets).*2,'ko');
hold on
plot(maxTargetsTrialN,miss(maxTargets),'mo');
plot(maxTargetsTrialN,fa(maxTargets).*0.5,'co');
figXAxis([],'Trial Number',[0 length(fa)])
figYAxis([],'True/False',[-1 3])
figAxForm([],0)
title({sprintf('%s Trials',num2str(length(fa)))})

if ~isempty(trialRange)
    visTargets = visTargets(trialRange);
    audTargets = audTargets(trialRange);
    invTrialInd = invTrialInd(trialRange);
    invVisTargets = invVisTargets(trialRange);
    invAudTargets = invAudTargets(trialRange);
    valOutcome = valOutcome(trialRange);
    invOutcome = invOutcome(trialRange);
end

hit = strcmp(valOutcome,'success');
miss = strcmp(valOutcome,'ignore');
fa = strcmp(valOutcome,'failure');
pctEarly = round(sum(fa)./length(valOutcome)*100);

invHit = strcmp(invOutcome,'FA');
invMiss = strcmp(invOutcome,'CR');
if length(invHit) ~= length(hit)
    extraTrials = length(hit) - length(invHit);
    invHit = cat(2,invHit,false(1,extraTrials));
    invMiss = cat(2,invMiss,false(1,extraTrials));
end

visTrials = visTargets > 0;
audTrials = audTargets > 0;
invVisTrials = invVisTargets > 0;
invAudTrials = invAudTargets > 0;

oris = unique(visTargets(visTrials));
invOris = unique(invVisTargets(invVisTrials));
amps = unique(audTargets(audTrials));
invAmps = unique(invAudTargets(invAudTrials));

nVisStim = length(oris);
nAudStim = length(amps);

nVisHit = zeros(1,nVisStim);
nVisMiss = zeros(1,nVisStim);
nInvVisHit = zeros(1,nVisStim);
nInvVisMiss = zeros(1,nVisStim);
nAudHit = zeros(1,nAudStim);
nAudMiss = zeros(1,nAudStim);
nInvAudHit = zeros(1,nAudStim);
nInvAudMiss = zeros(1,nAudStim);
for i = 1:nVisStim
    nVisHit(i) = sum(visTargets == oris(i) & hit);
    nVisMiss(i) = sum(visTargets == oris(i) & miss);
    if any(ismember(invOris,oris(i)))
        ind = invOris == oris(i);
        nInvVisHit(i) = sum(invVisTargets == invOris(ind) & invHit);
        nInvVisMiss(i) = sum(invVisTargets == invOris(ind) & invMiss);
    end    
end
for i = 1:nAudStim
    nAudHit(i) = sum(audTargets == amps(i) & hit);
    nAudMiss(i) = sum(audTargets == amps(i) & miss);  
    if any(ismember(invAmps,amps(i)))
        ind = invAmps == amps(i);
        nInvAudHit(i) = sum(invAudTargets == invAmps(ind) & invHit);
        nInvAudMiss(i) = sum(invAudTargets == invAmps(ind) & invMiss);
    end
end

pctInvVis = round((sum(nInvVisHit+nInvVisMiss)./...
    sum(nInvVisHit+nInvVisMiss+nVisHit+nVisMiss))*100);
pctInvAud = round((sum(nInvAudHit+nInvAudMiss)./...
    sum(nInvAudHit+nInvAudMiss+nAudHit+nAudMiss))*100);

invVisTargetInd = ismember(oris,invOris);
invAudTargetInd = ismember(amps,invAmps);

visHR = nVisHit./(nVisHit+nVisMiss);
invVisHR = nInvVisHit./(nInvVisHit+nInvVisMiss);
audHR = nAudHit./(nAudHit+nAudMiss);
invAudHR = nInvAudHit./(nInvAudHit+nInvAudMiss);

visInd = ~isnan(visHR);
if sum(visInd) >= 3
    visFit = weibullFitLG(oris(visInd), visHR(visInd),1, 0, {'nTrials',...
        nVisHit(visInd)+nVisMiss(visInd)});
else
    visFit = [];
end

audInd = ~isnan(audHR);
if sum(audInd) >= 3
    audFit = weibullFitLG(amps(audInd), audHR(audInd),1, 0, {'nTrials',...
        nAudHit(audInd)+nAudMiss(audInd)});
else
    audFit = [];
end

maxI = max(oris);
minI = min(oris);
visXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
visLevels_lim = [minI-10 maxI+10];
visLevels_label = [10 100];
maxI = max(amps);
minI = min(amps);
audXGrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
audLevels_lim = [minI/2 maxI+0.1];
audLevels_label = [0.001 0.01 0.1];

HR_lim = [0 110];
HR_label = 0:20:100;

suptitle({sprintf('%s-%s',ms,dt);sprintf('%s pct early',num2str(pctEarly))})
subplot 223
x = oris;
y = visHR*100;
h = plot(x,y,'ko');
hold on
x = invOris;
y = invVisHR(invVisTargetInd)*100;
h = plot(x,y,'co-');
if sum(visInd >= 3)
    h = line(visXGrid,visFit.modelFun(visFit.coefEsts,visXGrid)*100);
    h.Color = 'k';
    vline(visFit.thresh,'k--',sprintf('threshold = %s',num2str(visFit.thresh)));
end
ax = gca;
ax.XScale = 'log';
figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title({'Visual Trials';sprintf('%s pct inv',num2str(pctInvVis))})

subplot 224
x = amps;
y = audHR*100;
h = plot(x,y,'ko');
hold on
x = invAmps;
y = invAudHR(invAudTargetInd)*100;
h = plot(x,y,'co-');
if sum(audInd >= 3)
    h = line(audXGrid,audFit.modelFun(audFit.coefEsts,audXGrid)*100);
    h.Color = 'k';
    vline(audFit.thresh,'k--',sprintf('threshold = %s',num2str(audFit.thresh)));
end
ax = gca;
ax.XScale = 'log';
figXAxis([],'Amplitude (% Max)',audLevels_lim,audLevels_label,audLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title({'Auditory Trials';sprintf('%s pct inv',num2str(pctInvAud))})

if doSave
    print(fullfile(fnout,sprintf('%s-%s-HR',ms,dt)),'-dpdf','-fillpage')
end
%%
end
% end