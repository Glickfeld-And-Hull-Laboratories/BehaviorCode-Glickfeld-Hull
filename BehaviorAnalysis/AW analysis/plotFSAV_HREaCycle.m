% prevTrialAud = 1;
% prevTrialVis = 0;

valVisTrials = visTargets_all > 1;
valAudTrials = audTargets_all > 0;
invVisTrials = invVisTargets_all > 1;
invAudTrials = invAudTargets_all > 0;

cycles = unique(nCycles_all);
ncyc = length(cycles);

hitsPerCycle_vis = nan(1,ncyc);
missesPerCycle_vis = nan(1,ncyc);
invHitsPerCycle_vis = nan(1,ncyc);
invMissesPerCycle_vis = nan(1,ncyc);
hitsPerCycle_aud = nan(1,ncyc);
missesPerCycle_aud = nan(1,ncyc);
invHitsPerCycle_aud = nan(1,ncyc);
invMissesPerCycle_aud = nan(1,ncyc);

hitsPerCycle_prevVis = nan(1,ncyc);
missesPerCycle_prevVis = nan(1,ncyc);
hitsPerCycle_prevAud = nan(1,ncyc);
missesPerCycle_prevAud = nan(1,ncyc);

invHitsPerCycle_prevVis = nan(1,ncyc);
invMissesPerCycle_prevVis = nan(1,ncyc);
invHitsPerCycle_prevAud = nan(1,ncyc);
invMissesPerCycle_prevAud = nan(1,ncyc);
for icyc = 1:ncyc
    hitsPerCycle_vis(icyc) = sum(successIx_all & valVisTrials & nCycles_all == icyc);
    missesPerCycle_vis(icyc) = sum(missedIx_all & valVisTrials & nCycles_all == icyc);
    invHitsPerCycle_vis(icyc) = sum(invHitIx_all & invVisTrials & nCycles_all == icyc);
    invMissesPerCycle_vis(icyc) = sum(invMissIx_all & invVisTrials & nCycles_all == icyc);    
    
    hitsPerCycle_aud(icyc) = sum(successIx_all & valAudTrials & nCycles_all == icyc);
    missesPerCycle_aud(icyc) = sum(missedIx_all & valAudTrials & nCycles_all == icyc);
    invHitsPerCycle_aud(icyc) = sum(invHitIx_all & invAudTrials & nCycles_all == icyc);
    invMissesPerCycle_aud(icyc) = sum(invMissIx_all & invAudTrials & nCycles_all == icyc);    
    
    hitsPerCycle_prevVis(icyc) = sum(...
        successIx_all & valVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialVis);
    missesPerCycle_prevVis(icyc) = sum(...
        missedIx_all & valVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialVis);
    hitsPerCycle_prevAud(icyc) = sum(...
        successIx_all & valVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialAud);
    missesPerCycle_prevAud(icyc) = sum(...
        missedIx_all & valVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialAud);

    invHitsPerCycle_prevVis(icyc) = sum(...
        invHitIx_all & invVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialVis);
    invMissesPerCycle_prevVis(icyc) = sum(...
        invMissIx_all & invVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialVis);
    invHitsPerCycle_prevAud(icyc) = sum(...
        invHitIx_all & invVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialAud);
    invMissesPerCycle_prevAud(icyc) = sum(...
        invMissIx_all & invVisTrials & nCycles_all == icyc & prevTrial_all == prevTrialAud);
end

nValTrialsPerCycle_vis = hitsPerCycle_vis+missesPerCycle_vis;
[valHRPerCycle_vis, val95PerCycle_vis] = binofit(hitsPerCycle_vis,...
    nValTrialsPerCycle_vis);
nInvTrialsPerCycle_vis = invHitsPerCycle_vis+invMissesPerCycle_vis;
[invHRPerCycle_vis, inv95PerCycle_vis] = binofit(invHitsPerCycle_vis,...
    nInvTrialsPerCycle_vis);

nValTrialsPerCycle_aud = hitsPerCycle_aud+missesPerCycle_aud;
[valHRPerCycle_aud, val95PerCycle_aud] = binofit(hitsPerCycle_aud,...
    nValTrialsPerCycle_aud);
nInvTrialsPerCycle_aud = invHitsPerCycle_aud+invMissesPerCycle_aud;
[invHRPerCycle_aud, inv95PerCycle_aud] = binofit(invHitsPerCycle_aud,...
    nInvTrialsPerCycle_aud);

cyclesPlotInd_vis = nValTrialsPerCycle_vis >= minTrN_all &...
    nInvTrialsPerCycle_vis >= minTrN_all;
cyclesPlotInd_aud = nValTrialsPerCycle_aud >= minTrN_all &...
    nInvTrialsPerCycle_aud >= minTrN_all;
% cyclesPlotInd_aud = cyclesPlotInd_vis;

cyclesPlot_vis = cycles(cyclesPlotInd_vis);
cyclesPlot_aud = cycles(cyclesPlotInd_aud);

figure
subplot 121
h = errorbar(cyclesPlot_vis,valHRPerCycle_vis(cyclesPlotInd_vis),...
    valHRPerCycle_vis(cyclesPlotInd_vis) - val95PerCycle_vis(cyclesPlotInd_vis,1)',...
    val95PerCycle_vis(cyclesPlotInd_vis,2)' - valHRPerCycle_vis(cyclesPlotInd_vis),'ko');
h.MarkerFaceColor = [1 1 1];
hold on
h = errorbar(cyclesPlot_vis,invHRPerCycle_vis(cyclesPlotInd_vis),...
    invHRPerCycle_vis(cyclesPlotInd_vis) - inv95PerCycle_vis(cyclesPlotInd_vis,1)',...
    inv95PerCycle_vis(cyclesPlotInd_vis,2)' - invHRPerCycle_vis(cyclesPlotInd_vis),'co');
h.MarkerFaceColor = 'c';
figXAxis([],'Stimulus Number',[cyclesPlot_vis(1)-1 cyclesPlot_vis(end)+1],...
    cyclesPlot_vis,cyclesPlot_vis)
figYAxis([],'Hit Rate',[0 1])
figAxForm([])
title('visual trials')
subplot 122
h = errorbar(cyclesPlot_aud,valHRPerCycle_aud(cyclesPlotInd_aud),...
    valHRPerCycle_aud(cyclesPlotInd_aud) - val95PerCycle_aud(cyclesPlotInd_aud,1)',...
    val95PerCycle_aud(cyclesPlotInd_aud,2)' - valHRPerCycle_aud(cyclesPlotInd_aud),'ko');
h.MarkerFaceColor = [1 1 1];
hold on
h = errorbar(cyclesPlot_aud,invHRPerCycle_aud(cyclesPlotInd_aud),...
    invHRPerCycle_aud(cyclesPlotInd_aud) - inv95PerCycle_aud(cyclesPlotInd_aud,1)',...
    inv95PerCycle_aud(cyclesPlotInd_aud,2)' - invHRPerCycle_aud(cyclesPlotInd_aud),'co');
h.MarkerFaceColor = 'c';
figXAxis([],'Stimulus Number',[cyclesPlot_aud(1)-1 cyclesPlot_aud(end)+1],...
    cyclesPlot_aud,cyclesPlot_aud)
figYAxis([],'Hit Rate',[0 1])
figAxForm([])
title('auditory trials')
print(fullfile(fnout,'allTrainExpt_HReaCycle'),'-dpdf')

% %%
% minTrN_prev = 10;
% 
% nValTrialsPerCycle_prevVis = hitsPerCycle_prevVis+missesPerCycle_prevVis;
% [valHRPerCycle_prevVis, val95PerCycle_prevVis] = binofit(hitsPerCycle_prevVis,...
%     nValTrialsPerCycle_prevVis);
% nValTrialsPerCycle_prevAud = hitsPerCycle_prevAud+missesPerCycle_prevAud;
% [valHRPerCycle_prevAud, val95PerCycle_prevAud] = binofit(hitsPerCycle_prevAud,...
%     nValTrialsPerCycle_prevAud);
% 
% nInvTrialsPerCycle_prevVis = invHitsPerCycle_prevVis+invMissesPerCycle_prevVis;
% [invHRPerCycle_prevVis, inv95PerCycle_prevVis] = binofit(invHitsPerCycle_prevVis,...
%     nInvTrialsPerCycle_prevVis);
% nInvTrialsPerCycle_prevAud = invHitsPerCycle_prevAud+invMissesPerCycle_prevAud;
% [invHRPerCycle_prevAud, inv95PerCycle_prevAud] = binofit(invHitsPerCycle_prevAud,...
%     nInvTrialsPerCycle_prevAud);
% 
% cyclesPlotInd = nValTrialsPerCycle_prevVis >=minTrN_prev & ...
%     nValTrialsPerCycle_prevAud >=minTrN_prev & nInvTrialsPerCycle_prevVis...
%     >=minTrN_prev & nInvTrialsPerCycle_prevAud;
% 
% cyclesPlot_prev = cycles(cyclesPlotInd);
% 
% figure
% subplot 121
% h = errorbar(cyclesPlot_prev,valHRPerCycle_prevVis(cyclesPlotInd),...
%     valHRPerCycle_prevVis(cyclesPlotInd) - val95PerCycle_prevVis(cyclesPlotInd,1)',...
%     val95PerCycle_prevVis(cyclesPlotInd,2)' - valHRPerCycle_prevVis(cyclesPlotInd),'ko');
% h.MarkerFaceColor = [1 1 1];
% hold on
% 
% h = errorbar(cyclesPlot_prev,valHRPerCycle_prevAud(cyclesPlotInd),...
%     valHRPerCycle_prevAud(cyclesPlotInd) - val95PerCycle_prevAud(cyclesPlotInd,1)',...
%     val95PerCycle_prevAud(cyclesPlotInd,2)' - valHRPerCycle_prevAud(cyclesPlotInd),'ko');
% h.MarkerFaceColor = [1 1 1];
% h.Color = [0.5 0.5 0.5];
% 
% figXAxis([],'Stimulus Number',[cyclesPlot_prev(1)-1 cyclesPlot_prev(end)+1],...
%     cyclesPlot_prev,cyclesPlot_prev)
% figYAxis([],'Hit Rate',[0 1])
% figAxForm([])
% title('valid visual trials')
% legend({'prev vis';'prev aud'},'location','southwest')
% 
% subplot 122
% h = errorbar(cyclesPlot_prev,invHRPerCycle_prevVis(cyclesPlotInd),...
%     invHRPerCycle_prevVis(cyclesPlotInd) - inv95PerCycle_prevVis(cyclesPlotInd,1)',...
%     inv95PerCycle_prevVis(cyclesPlotInd,2)' - invHRPerCycle_prevVis(cyclesPlotInd),'co');
% h.MarkerFaceColor = [1 1 1];
% h.Color = [0 0.5 0.5];
% hold on
% 
% h = errorbar(cyclesPlot_prev,invHRPerCycle_prevAud(cyclesPlotInd),...
%     invHRPerCycle_prevAud(cyclesPlotInd) - inv95PerCycle_prevAud(cyclesPlotInd,1)',...
%     inv95PerCycle_prevAud(cyclesPlotInd,2)' - invHRPerCycle_prevAud(cyclesPlotInd),'co');
% h.MarkerFaceColor = [1 1 1];
% 
% figXAxis([],'Stimulus Number',[cyclesPlot_prev(1)-1 cyclesPlot_prev(end)+1],...
%     cyclesPlot_prev,cyclesPlot_prev)
% figYAxis([],'Hit Rate',[0 1])
% figAxForm([])
% title('invalid visual trials')
% legend({'prev vis';'prev aud'},'location','southwest')
% 
% print(fullfile(fnout,'allTrainExpt_HReaCyclePrevTrial'),'-dpdf')