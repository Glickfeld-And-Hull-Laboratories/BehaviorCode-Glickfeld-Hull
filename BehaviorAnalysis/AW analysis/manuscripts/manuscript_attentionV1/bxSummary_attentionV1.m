clear all
close all
rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
datasetStr = 'FSAV_V1_100ms';
fnout = fullfile(rc.caOutputDir,datasetStr,'behavior');
if ~exist(fnout,'dir')
    mkdir(fnout)
end
bxImagingDatasets;
imagingDatasets = expt;
clear expt

rewardTrainedMice = {'613';'614'};
%%
minTrN_expt = 2;
minTrN_ms = 5;
minTrN_all = 5;
thresh = 30;
invVisHRCutoff = 23;
invAudHRCutoff = 0.01;
early_cutoff = 0.5;
lapse_cutoff = 0.85;
thresh_ms = 2000;
nTimes = 8;
nBins = 6;
visualTrials = 1;
auditoryTrials = 2;
val = 1;
inv = 2;
oriBin = 1;
timeBin = 2;

nexp = xd.nRows;
mice = unique(xd.Subject);
nmice = length(mice);
colors_mice = parula(nmice);
nexp_img = size(imagingDatasets,2);
% expCatch = logical(cell2mat({expt.catch}));


doExptStruct = 0;
if doExptStruct
    bxExptStruct_FSAV_train
    bxExptStruct_FSAV_img
else
    load(fullfile(fnout,'bxExpMat_training.mat'))
    load(fullfile(fnout,'bxExpMat_imaging.mat'))
end

amp_edges = [0.003 0.007 0.02 0.05 0.14 0.37 1];
% amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];


%% analyze behavior data from each experiment
ms = struct;
for im = 1:nmice
    ms(im).sn = mice(im);
end
visTargets_all = [];
invVisTargets_all = [];
invAudTargets_all = [];
invCycles_all = [];
audTargets_all = [];
successIx_all = [];
missedIx_all = [];
invHitIx_all = [];
invMissIx_all = [];
nCycles_all = [];
trLengthMs_all = [];
prevTrial_all = [];
prev2Trial_all = [];
prevVisTarget_all = [];

pctInv_expt = zeros(1,nexp + nexp_img);
invMatchHR_expt = zeros(1,nexp + nexp_img);
valMatchHR_expt = zeros(1,nexp + nexp_img);
doPlot = 0;
ntr1 = 0;
ntr2 = 0;
expAnalyzed = ones(1,nexp_img);
for iexp = 1:nexp
    if isnan(bxExp(iexp).invType)
        continue
    else
        ntr1 = ntr1 + length(bxExp(iexp).sIx);
%         if early_mat(iexp) > early_cutoff | HR_ori_mat(iexp) < lapse_cutoff | HR_amp_mat(iexp) < lapse_cutoff | stimOnTime100 == 0
%             continue
%         end
        visHRCutoffPass = visHitRate{iexp}(end-1:end) < lapse_cutoff;
        if length(audHitRate{iexp}) < 2
            audHRCutoffPass = audHitRate{iexp};
        else
            audHRCutoffPass = audHitRate{iexp}(end-1:end) < lapse_cutoff;
        end
        if early_mat(iexp) > early_cutoff | ...
            all(visHRCutoffPass) | all(audHRCutoffPass)
            expAnalyzed(iexp) = 0;
            continue
        end
        ntr2 = ntr2 + length(bxExp(iexp).sIx);
        sn = bxExp(iexp).sn;
        sn_ind = find(mice == sn);
        cycLengthMs = double(bxExp(iexp).tOn + bxExp(iexp).tOff);

        visTargets = chop(double(bxExp(iexp).tVisTargets),2);
        audTargets = chop(double(bxExp(iexp).tAudTargets),2);

        nt = length(visTargets);
        
        trType = visTargets > 0;
        trType_shift = [NaN trType];
        prevTrial = trType_shift(1:length(trType));
        trType_shift = [NaN NaN trType];
        prev2Trial = trType_shift(1:length(trType));

        visTargets_shift = [NaN visTargets];
        prevVisTarget = visTargets_shift(1:length(visTargets));

        successIx = double(bxExp(iexp).sIx);
        missedIx = double(bxExp(iexp).mIx);
        invHitIx = double(bxExp(iexp).invHitIx);
        invMissIx = double(bxExp(iexp).invMissIx);
        nCyc = double(bxExp(iexp).trLength);
        trLengthMs = nCyc.*cycLengthMs;

        % sort vis catch and aud catch days, match invalid to valid target
        % types
        try
            if strcmp(bxExp(iexp).invType,'vis')
                invVisTargets = chop(double(bxExp(iexp).tInvTargets),2);
                invAudTargets = nan(1,nt);
                invTrials = invVisTargets > 0;
                invTars = round(unique(invVisTargets(invTrials)));
                matchValTrials = matchTrialsBx(invVisTargets,invTars,invTrials,visTargets);
                pctInv_expt(iexp) = sum(invTrials)./nt;
                invMatchHR_expt(iexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
                valMatchHR_expt(iexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));
                
                invHR_expt = nan(1:length(invTars));
                valHR_expt = nan(1:length(invTars));
                for i = 1:length(invTars)
                    invHR_expt(i) = sum(invHitIx & invVisTargets == invTars(i))./ ...
                        sum((invHitIx | invMissIx) & invVisTargets == invTars(i));
                    valHR_expt(i) = sum(successIx & visTargets == invTars(i))./...
                        sum((successIx | missedIx) & visTargets == invTars(i));
                end
                
                nInvTrials(iexp) = sum(invTrials);
                invCycles = bxExp(iexp).invTrLength;
            elseif strcmp(bxExp(iexp).invType,'aud')
                invVisTargets = nan(1,nt);
                invAudTargets = chop(double(bxExp(iexp).tInvTargets),2);  
                invTrials = invAudTargets > 0;
                invTars = unique(invAudTargets(invTrials));  
                matchValTrials = matchTrialsBx(invAudTargets,invTars,invTrials,audTargets); 
                pctInv_expt(iexp) = sum(invTrials)./nt;
                invMatchHR_expt(iexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
                valMatchHR_expt(iexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));  
                
                invHR_expt = nan(1:length(invTars));
                valHR_expt = nan(1:length(invTars));
                for i = 1:length(invTars)
                    invHR_expt(i) = sum(invHitIx & invAudTargets == invTars(i))./ ...
                        sum((invHitIx | invMissIx) & invAudTargets == invTars(i));
                    valHR_expt(i) = sum(successIx & audTargets == invTars(i))./...
                        sum((successIx | missedIx) & audTargets == invTars(i));
                end
                
                nInvTrials(iexp) = sum(invTrials);
                invCycles = bxExp(iexp).invTrLength;
            else
                invVisTargets = nan(1,nt);
                invAudTargets = nan(1,nt);
                invCycles = nan(1,nt);
                pctInv_expt(iexp) = nan;
                invMatchHR_expt(iexp) = nan;
                valMatchHR_expt(iexp) = nan;
                nInvTrials(iexp) = nan;
                invHR_expt = nan;
                valHR_expt = nan;
                invTars = nan;
            end
        catch
            disp(iexp)
            invVisTargets = nan(1,nt);
            invAudTargets = nan(1,nt);
            invCycles = nan(1,nt);
            pctInv_expt(iexp) = nan;
            invMatchHR_expt(iexp) = nan;
            valMatchHR_expt(iexp) = nan;
            nInvTrials(iexp) = nan; 
            invHR_expt = nan;
            valHR_expt= nan;  
            invTars = nan
        end
        
        
        if strcmp(bxExp(iexp).invType,'vis')
            matches = getMatchedValidTrialIndex(bxExp(iexp).tVisTargets,...
                bxExp(iexp).tInvTargets);
        elseif strcmp(bxExp(iexp).invType,'aud')
            matches = getMatchedValidTrialIndex(bxExp(iexp).tAudTargets,...
                bxExp(iexp).tInvTargets);
        else
            matches = {nan};
        end
        matchedValidSamples = cell2mat(matches);
        matchedValidHit = successIx(matchedValidSamples);
        matchedValidMiss = missedIx(matchedValidSamples);
        matchedValidHitRate = sum(matchedValidHit)./sum(matchedValidHit | ...
            matchedValidMiss);
        matchedInvalidHitRate = sum(invHitIx)./sum(invHitIx | invMissIx);
        
        visTargets_all = cat(2,visTargets_all,visTargets);
        invVisTargets_all = cat(2,invVisTargets_all,invVisTargets);
        invAudTargets_all = cat(2,invAudTargets_all,invAudTargets);
        invCycles_all = cat(2,invCycles_all,invCycles);
        audTargets_all = cat(2,audTargets_all,audTargets);
        successIx_all = cat(2,successIx_all,successIx);
        missedIx_all = cat(2,missedIx_all,missedIx);
        invHitIx_all = cat(2,invHitIx_all,invHitIx);
        invMissIx_all = cat(2,invMissIx_all,invMissIx);
        trLengthMs_all = cat(2,trLengthMs_all,trLengthMs);
        nCycles_all = cat(2,nCycles_all,nCyc);
        prevTrial_all = cat(2,prevTrial_all,prevTrial);
        prev2Trial_all = cat(2,prev2Trial_all,prev2Trial);
        prevVisTarget_all = cat(2,prevVisTarget_all,prevVisTarget);

        if isfield(ms,'visTargets')           
            ms(sn_ind).visTargets = cat(2,ms(sn_ind).visTargets,visTargets);
            ms(sn_ind).invVisTargets = cat(2,ms(sn_ind).invVisTargets,invVisTargets);
            ms(sn_ind).invAudTargets = cat(2,ms(sn_ind).invAudTargets,invAudTargets);
            ms(sn_ind).invCycles = cat(2,ms(sn_ind).invCycles,invCycles);
            ms(sn_ind).audTargets = cat(2,ms(sn_ind).audTargets,audTargets);
            ms(sn_ind).successIx = cat(2,ms(sn_ind).successIx,successIx);
            ms(sn_ind).missedIx = cat(2,ms(sn_ind).missedIx,missedIx);
            ms(sn_ind).invHitIx = cat(2,ms(sn_ind).invHitIx,invHitIx);
            ms(sn_ind).invMissIx = cat(2,ms(sn_ind).invMissIx,invMissIx);
            ms(sn_ind).nCycles = cat(2,ms(sn_ind).nCycles,nCyc);
            ms(sn_ind).trLengthMs = cat(2,ms(sn_ind).trLengthMs,trLengthMs);
            ms(sn_ind).invType = cat(2,ms(sn_ind).invType,{bxExp(iexp).invType});
            ms(sn_ind).valMatchedHR = cat(2, ms(sn_ind).valMatchedHR,...
                matchedValidHitRate);
            ms(sn_ind).invMatchedHR = cat(2, ms(sn_ind).invMatchedHR,...
                matchedInvalidHitRate);
            ms(sn_ind).valHRExpt = cat(2, ms(sn_ind).valHRExpt,{valHR_expt});
            ms(sn_ind).invHRExpt = cat(2, ms(sn_ind).invHRExpt,{invHR_expt});
            ms(sn_ind).invTargetsExpt = cat(2,ms(sn_ind).invTargetsExpt,{invTars});
            
        else
            ms(sn_ind).visTargets = visTargets;
            ms(sn_ind).invVisTargets = invVisTargets;
            ms(sn_ind).invAudTargets = invAudTargets;
            ms(sn_ind).invCycles = invCycles;
            ms(sn_ind).audTargets = audTargets;
            ms(sn_ind).successIx = successIx;
            ms(sn_ind).missedIx = missedIx;
            ms(sn_ind).invHitIx = invHitIx;
            ms(sn_ind).invMissIx = invMissIx;  
            ms(sn_ind).nCycles = nCyc;
            ms(sn_ind).trLengthMs = trLengthMs;
            ms(sn_ind).invType = {bxExp(iexp).invType};
            ms(sn_ind).valMatchedHR = matchedValidHitRate;
            ms(sn_ind).invMatchedHR = matchedInvalidHitRate;
            ms(sn_ind).valHRExpt = {valHR_expt};
            ms(sn_ind).invHRExpt = {invHR_expt};
            ms(sn_ind).invTargetsExpt = {invTars};
        end

        if length(ms(sn_ind).visTargets) ~= length(ms(sn_ind).invHitIx)
            error(['error on expt ' num2str(iexp)])
        end
        
%         bxSumExptHR_FSAV_img
    end
end

expAnalyzed_img = ones(1,nexp_img);
for iexp = 1:nexp_img    
    if isnan(bxImgExp(iexp).invType)
        continue
    else
%         if early_mat_img(iexp) > early_cutoff | HR_ori_mat_img(iexp) < lapse_cutoff | HR_amp_mat_img(iexp) < lapse_cutoff | stimOnTime100_img(iexp) == 0
%             continue
%         end
        visHRCutoffPass_img = visHitRate_img{iexp}(end-1:end) < lapse_cutoff;
        audHRCutoffPass_img = audHitRate_img{iexp}(end-1:end) < lapse_cutoff;
        if early_mat_img(iexp) > early_cutoff | ...
            all(visHRCutoffPass_img) | all(audHRCutoffPass_img)
            expAnalyzed(iexp) = 0;
            continue
        end
        sn = bxImgExp(iexp).sn;
        sn_ind = find(mice == str2num(sn));
        cycLengthMs = double(bxImgExp(iexp).tOn + bxImgExp(iexp).tOff);
        visTargets = chop(double(bxImgExp(iexp).tVisTargets),2);
        audTargets = chop(double(bxImgExp(iexp).tAudTargets),2);

        nt = length(visTargets);
        
        trType = visTargets > 0;
        trType_shift = [NaN trType];
        prevTrial = trType_shift(1:length(trType));
        trType_shift = [NaN NaN trType];
        prev2Trial = trType_shift(1:length(trType));

        visTargets_shift = [NaN visTargets];
        prevVisTarget = visTargets_shift(1:length(visTargets));

        successIx = double(bxImgExp(iexp).sIx);
        missedIx = double(bxImgExp(iexp).mIx);
        invHitIx = double(bxImgExp(iexp).invHitIx);
        invMissIx = double(bxImgExp(iexp).invMissIx);
        nCyc = double(bxImgExp(iexp).trLength);
        trLengthMs = nCyc.*cycLengthMs;

        % sort vis catch and aud catch days, match invalid to valid target
        % types
        try
            if strcmp(bxImgExp(iexp).invType,'vis')
                invVisTargets = chop(double(bxImgExp(iexp).tInvTargets),2);
                invAudTargets = nan(1,nt);
                invTrials = invVisTargets > 0;
                invTars = round(unique(invVisTargets(invTrials)));
                matchValTrials = matchTrialsBx(invVisTargets,invTars,invTrials,visTargets);
                pctInv_expt(iexp+nexp) = sum(invTrials)./nt;
                invMatchHR_expt(iexp+nexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
                valMatchHR_expt(iexp+nexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));
                nInvTrials(iexp+nexp) = sum(invTrials);
                invCycles = bxImgExp(iexp).invTrLength;
                
                invHR_expt = nan(1:length(invTars));
                valHR_expt = nan(1:length(invTars));
                for i = 1:length(invTars)
                    invHR_expt(i) = sum(invHitIx & invVisTargets == invTars(i))./ ...
                        sum((invHitIx | invMissIx) & invVisTargets == invTars(i));
                    valHR_expt(i) = sum(successIx & visTargets == invTars(i))./...
                        sum((successIx | missedIx) & visTargets == invTars(i));
                end
            elseif strcmp(bxImgExp(iexp).invType,'aud')
                invVisTargets = nan(1,nt);
                invAudTargets = chop(double(bxImgExp(iexp).tInvTargets),2);  
                invTrials = invAudTargets > 0;
                invTars = unique(invAudTargets(invTrials));  
                matchValTrials = matchTrialsBx(invAudTargets,invTars,invTrials,audTargets); 
                pctInv_expt(iexp+nexp) = sum(invTrials)./nt;
                invMatchHR_expt(iexp+nexp) = sum(successIx(invTrials))./(sum(invHitIx(invTrials)) + sum(invMissIx(invTrials)));
                valMatchHR_expt(iexp+nexp) = sum(successIx(matchValTrials))./(sum(successIx(matchValTrials)) + sum(missedIx(matchValTrials)));  
                nInvTrials(iexp+nexp) = sum(invTrials);
                invCycles = bxImgExp(iexp).invTrLength;
                invHR_expt= nan(1:length(invTars));
                valHR_expt = nan(1:length(invTars));
                for i = 1:length(invTars)
                    invHR_expt(i) = sum(invHitIx & invAudTargets == invTars(i))./ ...
                        sum((invHitIx | invMissIx) & invAudTargets == invTars(i));
                    valHR_expt(i) = sum(successIx & audTargets == invTars(i))./...
                        sum((successIx | missedIx) & audTargets == invTars(i));
                end
            else
                invVisTargets = nan(1,nt);
                invAudTargets = nan(1,nt);
                invCycles = nan(1,nt);
                pctInv_expt(iexp+nexp) = nan;
                invMatchHR_expt(iexp+nexp) = nan;
                valMatchHR_expt(iexp+nexp) = nan;
                nInvTrials(iexp+nexp) = nan;
                invHR_expt = nan;
                valHR_expt = nan;
            end
        catch
            invVisTargets = nan(1,nt);
            invAudTargets = nan(1,nt);
            invCycles = nan(1,nt);
            pctInv_expt(iexp+nexp) = nan;
            invMatchHR_expt(iexp+nexp) = nan;
            valMatchHR_expt(iexp+nexp) = nan;
            nInvTrials(iexp+nexp) = nan;  
            invHR_expt = nan;
            valHR_expt = nan;          
        end
        
        if strcmp(bxImgExp(iexp).invType,'vis')
            matches = getMatchedValidTrialIndex(bxImgExp(iexp).tVisTargets,...
                bxImgExp(iexp).tInvTargets);
        elseif strcmp(bxImgExp(iexp).invType,'aud')
            matches = getMatchedValidTrialIndex(bxImgExp(iexp).tAudTargets,...
                bxImgExp(iexp).tInvTargets);
        else
            matches = {nan};
        end
        matchedValidSamples = cell2mat(matches);
        matchedValidHit = successIx(matchedValidSamples);
        matchedValidMiss = missedIx(matchedValidSamples);
        matchedValidHitRate = sum(matchedValidHit)./sum(matchedValidHit | ...
            matchedValidMiss);
        matchedInvalidHitRate = sum(invHitIx)./sum(invHitIx | invMissIx);
        
        visTargets_all = cat(2,visTargets_all,visTargets);
        invVisTargets_all = cat(2,invVisTargets_all,invVisTargets);
        invAudTargets_all = cat(2,invAudTargets_all,invAudTargets);
        invCycles_all = cat(2,invCycles_all,invCycles);
        audTargets_all = cat(2,audTargets_all,audTargets);
        successIx_all = cat(2,successIx_all,successIx);
        missedIx_all = cat(2,missedIx_all,missedIx);
        invHitIx_all = cat(2,invHitIx_all,invHitIx);
        invMissIx_all = cat(2,invMissIx_all,invMissIx);
        trLengthMs_all = cat(2,trLengthMs_all,trLengthMs);
        nCycles_all = cat(2,nCycles_all,nCyc);
        prevTrial_all = cat(2,prevTrial_all,prevTrial);
        prev2Trial_all = cat(2,prev2Trial_all,prev2Trial);
        prevVisTarget_all = cat(2,prevVisTarget_all,prevVisTarget);

        if isfield(ms,'visTargets')           
            ms(sn_ind).visTargets = cat(2,ms(sn_ind).visTargets,visTargets);
            ms(sn_ind).invVisTargets = cat(2,ms(sn_ind).invVisTargets,invVisTargets);
            ms(sn_ind).invAudTargets = cat(2,ms(sn_ind).invAudTargets,invAudTargets);
            ms(sn_ind).invCycles = cat(2,ms(sn_ind).invCycles,invCycles);
            ms(sn_ind).audTargets = cat(2,ms(sn_ind).audTargets,audTargets);
            ms(sn_ind).successIx = cat(2,ms(sn_ind).successIx,successIx);
            ms(sn_ind).missedIx = cat(2,ms(sn_ind).missedIx,missedIx);
            ms(sn_ind).invHitIx = cat(2,ms(sn_ind).invHitIx,invHitIx);
            ms(sn_ind).invMissIx = cat(2,ms(sn_ind).invMissIx,invMissIx);
            ms(sn_ind).nCycles = cat(2,ms(sn_ind).nCycles,nCyc);
            ms(sn_ind).trLengthMs = cat(2,ms(sn_ind).trLengthMs,trLengthMs);
            ms(sn_ind).invType = cat(2,ms(sn_ind).invType,{bxImgExp(iexp).invType});
            ms(sn_ind).valMatchedHR = cat(2, ms(sn_ind).valMatchedHR,...
                matchedValidHitRate);
            ms(sn_ind).invMatchedHR = cat(2, ms(sn_ind).invMatchedHR,...
                matchedInvalidHitRate);
            ms(sn_ind).valHRExpt = cat(2, ms(sn_ind).valHRExpt,{valHR_expt});
            ms(sn_ind).invHRExpt = cat(2, ms(sn_ind).invHRExpt,{invHR_expt});
            ms(sn_ind).invTargetsExpt = cat(2,ms(sn_ind).invTargetsExpt,{invTars});
        else
            ms(sn_ind).visTargets = visTargets;
            ms(sn_ind).invVisTargets = invVisTargets;
            ms(sn_ind).invAudTargets = invAudTargets;
            ms(sn_ind).invCycles = invCycles;
            ms(sn_ind).audTargets = audTargets;
            ms(sn_ind).successIx = successIx;
            ms(sn_ind).missedIx = missedIx;
            ms(sn_ind).invHitIx = invHitIx;
            ms(sn_ind).invMissIx = invMissIx; 
            ms(sn_ind).nCycles = nCyc;
            ms(sn_ind).trLengthMs = trLengthMs;
            ms(sn_ind).invType = {bxImgExp(iexp).invType};
            ms(sn_ind).valMatchedHR = matchedValidHitRate;
            ms(sn_ind).invMatchedHR = matchedInvalidHitRate;
            ms(sn_ind).valHRExpt = {valHR_expt};
            ms(sn_ind).invHRExpt = {invHR_expt};
            ms(sn_ind).invTargetsExpt = {invTars};
        end

        if length(ms(sn_ind).visTargets) ~= length(ms(sn_ind).invHitIx)
            error(['error on expt ' num2str(iexp)])
        end
    end
end

%% behavior summary for each mouse
msName = cellfun(@num2str,{ms.sn},'unif',0);
nVisSessionsEaMs = zeros(1,nmice);
nAudSessionsEaMs = zeros(1,nmice);
for im = 1:nmice
    nVisSessionsEaMs(im) = sum(cellfun(@(x,y,z,a,b)...
        x == mice(im) & strcmp(y,'vis') & z < early_cutoff...
        & a > lapse_cutoff & b > lapse_cutoff, {bxExp.sn},{bxExp.invType},...
        num2cell(early_mat)',num2cell(HR_ori_mat)',num2cell(HR_amp_mat)'))...
        + sum(cellfun(@(x,y,z,a,b)...
        strcmp(x,num2str(mice(im))) & strcmp(y,'vis') & z < early_cutoff...
        & a > lapse_cutoff & b > lapse_cutoff, {bxImgExp.sn},{bxImgExp.invType},...
        num2cell(early_mat_img)',num2cell(HR_ori_mat_img)',num2cell(HR_amp_mat_img)'));  
    nAudSessionsEaMs(im) = sum(cellfun(@(x,y,z,a,b)...
        x == mice(im) & strcmp(y,'aud') & z < early_cutoff...
        & a > lapse_cutoff & b > lapse_cutoff, {bxExp.sn},{bxExp.invType},...
        num2cell(early_mat)',num2cell(HR_ori_mat)',num2cell(HR_amp_mat)'))...
        + sum(cellfun(@(x,y,z,a,b)...
        strcmp(x,num2str(mice(im))) & strcmp(y,'aud') & z < early_cutoff...
        & a > lapse_cutoff & b > lapse_cutoff, {bxImgExp.sn},{bxImgExp.invType},...
        num2cell(early_mat_img)',num2cell(HR_ori_mat_img)',num2cell(HR_amp_mat_img)'));    
end

msSummary = getEaMouseBxSummary(ms,nBins,minTrN_ms);

nVisValid = zeros(1,nmice);
nVisInvalid = zeros(1,nmice);
nAudValid = zeros(1,nmice);
nAudInvalid = zeros(1,nmice);
for im = 1:nmice
    ind = ~isnan(msSummary(im).av(visualTrials).binning(oriBin).cue(...
        val).hitRate);
    nVisValid(im) = sum(msSummary(im).av(visualTrials).binning(oriBin).cue(...
        val).nTrials(ind));
    ind = ~isnan(msSummary(im).av(visualTrials).binning(oriBin).cue(...
        inv).hitRate);
    nVisInvalid(im) = sum(msSummary(im).av(visualTrials).binning(oriBin).cue(...
        inv).nTrials(ind));
    ind = ~isnan(msSummary(im).av(auditoryTrials).binning(oriBin).cue(...
        val).hitRate);
    nAudValid(im) = sum(msSummary(im).av(auditoryTrials).binning(oriBin).cue(...
        val).nTrials(ind));
    ind = ~isnan(msSummary(im).av(auditoryTrials).binning(oriBin).cue(...
        inv).hitRate);
    nAudInvalid(im) = sum(msSummary(im).av(auditoryTrials).binning(oriBin).cue(...
        inv).nTrials(ind));
end

bxSummaryTable = table(msName',nVisSessionsEaMs',nVisValid',nVisInvalid',...
    nAudSessionsEaMs',nAudValid',nAudInvalid');
bxSummaryTable.Properties.VariableNames = {'mouse','nVisSessions',...
    'nVisValidIncl','nVisInvalidIncl','nAudSessions','nAudValidIncl',...
    'nAudInvalidIncl'};

writetable(bxSummaryTable,fullfile(fnout,'bxSummaryTable.xls'),'WriteRowNames',true)

%% plotting params and targets
HR_lim = [0 110];
HR_label = [0:20:110];
cueColors = {'k';'c'};
cueNames = {'Valid';'Invalid'};

diffColors = {'b','k'};
visDiffNames = {sprintf('Target <= %s', num2str(invVisHRCutoff)),...
                sprintf('Target >= %s', num2str(invVisHRCutoff))};
audDiffNames = {sprintf('Target <= %s', num2str(invAudHRCutoff)),...
                sprintf('Target >= %s', num2str(invAudHRCutoff))};

visTargets = unique(visTargets_all);
visTargets = visTargets(2:end);
audTargets = unique(audTargets_all);
audTargets = audTargets(2:end);

visLevels_lim = [min(visTargets)-1 110];
visLevels_label = [10 100];
audLevels_lim = [min(audTargets)-(0.5*min(audTargets)) 1.1];
audLevels_label = [0 0.001 0.01 0.1 1];
%% individual mouse plots
for im = 1:nmice
    figure;
    suptitle(msName{im})
    
    subplot 221
    for icue = 1:2
        y = msSummary(im).av(visualTrials).binning(oriBin).cue(icue).hitRate .*100;
        yErrL = (y - msSummary(im).av(visualTrials).binning(oriBin).cue(icue).hitRateCI(:,1)'.*100);
        yErrU = (msSummary(im).av(visualTrials).binning(oriBin).cue(icue).hitRateCI(:,2)'.*100) - y;
        x = msSummary(im).av(visualTrials).binning(oriBin).cue(icue).binTarget;
        xErr = msSummary(im).av(visualTrials).binning(oriBin).cue(icue).binTargetSte;
        h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
        h.MarkerFaceColor = cueColors{icue};
        h.Color = cueColors{icue};
        hold on
    end
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Visual Trials')
    
    subplot 222
    for icue = 1:2
        y = msSummary(im).av(auditoryTrials).binning(oriBin).cue(icue).hitRate .*100;
        yErrL = (y - msSummary(im).av(auditoryTrials).binning(oriBin).cue(icue).hitRateCI(:,1)'.*100);
        yErrU = (msSummary(im).av(auditoryTrials).binning(oriBin).cue(icue).hitRateCI(:,2)'.*100)-y;
        x = msSummary(im).av(auditoryTrials).binning(oriBin).cue(icue).binTarget;
        xErr = msSummary(im).av(auditoryTrials).binning(oriBin).cue(icue).binTargetSte;
        h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
        h.MarkerFaceColor = [1 1 1];
        h.Color = cueColors{icue};
        hold on
    end
    ax = gca;
    ax.XScale = 'log';
    figXAxis([],'Tone Volume (% of max)',audLevels_lim,audLevels_label,audLevels_label)
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Auditory Trials')
    
    subplot 223
    invTypeInd = strcmp(ms(im).invType,'vis');
    nanInd = cellfun(@(x) all(isnan(x)),ms(im).invHRExpt);
    targets = cell2mat(ms(im).invTargetsExpt(invTypeInd & ~nanInd));
    valHR = cell2mat(ms(im).valHRExpt(invTypeInd & ~nanInd))*100;
    invHR = cell2mat(ms(im).invHRExpt(invTypeInd & ~nanInd))*100;
    targetInd = targets <= invVisHRCutoff;
    x = valHR(targetInd);
    y = invHR(targetInd);
    scatter(x,y,[diffColors{1} 'o']);
    hold on
    leg(1) = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),...
        ste(x,2),ste(x,2),[diffColors{1} 'o']);
    leg(1).MarkerFaceColor = diffColors{1};
    targetInd = targets >= invVisHRCutoff;
    x = valHR(targetInd);
    y = invHR(targetInd);
    scatter(x,y,[diffColors{2} 'o']);
    leg(2) = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),...
        ste(x,2),ste(x,2),[diffColors{2} 'o']);
    leg(2).MarkerFaceColor = diffColors{2};
    plot(HR_label,HR_label,'k--')
    figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Visual Expts')
    legend(leg,visDiffNames,'location','west')
    
    subplot 224
    invTypeInd = strcmp(ms(im).invType,'aud');
    nanInd = cellfun(@(x) all(isnan(x)),ms(im).invHRExpt);
    targets = cell2mat(ms(im).invTargetsExpt(invTypeInd & ~nanInd));
    valHR = cell2mat(ms(im).valHRExpt(invTypeInd & ~nanInd))*100;
    invHR = cell2mat(ms(im).invHRExpt(invTypeInd & ~nanInd))*100;
    targetInd = targets <= invAudHRCutoff;
    x = valHR(targetInd);
    y = invHR(targetInd);
    scatter(x,y,[diffColors{1} 'o']);
    hold on
    leg(1) = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),...
        ste(x,2),ste(x,2),[diffColors{1} 'o']);
    leg(1).MarkerFaceColor = diffColors{1};
    targetInd = targets >= invAudHRCutoff;
    x = valHR(targetInd);
    y = invHR(targetInd);
    scatter(x,y,[diffColors{2} 'o']);
    leg(2) = errorbar(nanmean(x),nanmean(y),ste(y,2),ste(y,2),...
        ste(x,2),ste(x,2),[diffColors{2} 'o']);
    leg(2).MarkerFaceColor = diffColors{2};
    plot(HR_label,HR_label,'k--')
    figXAxis([],'Valid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figYAxis([],'Invalid Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm
    title('Auditory Expts')
    legend(leg,audDiffNames,'location','west')
    
    print(fullfile(fnout, [msName{im} '_HRsummary']),'-dpdf','-fillpage')
end

figure
x = 1:2;
for im = 1:nmice
    subplot 121
    hold on
    y = [msSummary(im).av(visualTrials).binning(oriBin).cue(val).matchedHR ...
        msSummary(im).av(visualTrials).binning(oriBin).cue(inv).matchedHR].*100;
    h = plot(x,y,'ko-');
    h.Color = [0.5 0.5 0.5];
    h.MarkerFaceColor = [1 1 1];
    figXAxis([],'Cue Type',[0 3],x,cueNames);
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm([],0)
    title('Visual Trials')

    subplot 122
    hold on
    y = [msSummary(im).av(auditoryTrials).binning(oriBin).cue(val).matchedHR ...
        msSummary(im).av(auditoryTrials).binning(oriBin).cue(inv).matchedHR].*100;
    h = plot(x,y,'ko-');
    h.Color = [0.5 0.5 0.5];
    h.MarkerFaceColor = [1 1 1];
    figXAxis([],'Cue Type',[0 3],x,cueNames);
    figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
    figAxForm([],0)
    title('Auditory Trials')
end

%% threshold for each mouse
msFits = cell(1,nmice);
msXGrid = cell(1,nmice);
for im = 1:nmice
    orientations = msSummary(im).av(visualTrials).binning(oriBin).cue(val).binTarget;
    HR = msSummary(im).av(visualTrials).binning(oriBin).cue(val).hitRate;
    nTrials = msSummary(im).av(visualTrials).binning(oriBin).cue(val).nTrials;

    msFits{1,im} = weibullFitLG(orientations, HR,1, 0, {'nTrials', nTrials});

    maxI = max(orientations);
    minI = min(orientations);
    msXGrid{im} = logspace(log10(minI*0.1),log10(maxI*1.5),100);
end

rewMiceInd = ismember(msName,rewardTrainedMice);
msThresh = nan(1,length(rewMiceInd));
for im = 1:nmice
    msThresh(im) = msFits{im}.thresh;
end

%%
visInvCutoff = 24;
msMatchedHR = struct;
for im = 1:nmice
    invVisTargets = ms(im).invVisTargets;
    invVisTrials = invVisTargets > 0;
    hit = ms(im).invHitIx;
    miss = ms(im).invMissIx;
    ind1 = invVisTrials & (invVisTargets < visInvCutoff);
    ind2 = invVisTrials & (invVisTargets > visInvCutoff);
    invVisHR = [sum(hit(ind1))./sum(hit(ind1) | miss(ind1)) ...
        sum(hit(ind2))./sum(hit(ind2) | miss(ind2))];
    msMatchedHR(im).invVisHR = invVisHR;
    
    visTargets = ms(im).visTargets;
    invVisTargets(isnan(invVisTargets)) = 0;
    matchedTrials = cell2mat(getMatchedValidTrialIndex(visTargets,...
        invVisTargets));
    matchedVisTargets = visTargets(matchedTrials);
    matchHit = ms(im).successIx(matchedTrials);
    matchMiss = ms(im).missedIx(matchedTrials);
    ind1 = matchedVisTargets < visInvCutoff;
    ind2 = matchedVisTargets > visInvCutoff;
    visHR = [sum(matchHit(ind1))./sum(matchHit(ind1) | matchMiss(ind1)) ...
        sum(matchHit(ind2))./sum(matchHit(ind2) | matchMiss(ind2))];
    msMatchedHR(im).visHR = visHR;
    msMatchedHR(im).sub = visHR - invVisHR;
end


figure;
subplot 221
x = ones(1,sum(rewMiceInd));
y = msThresh(rewMiceInd);
h = plot(x,y,'ko');
h.MarkerFaceColor = 'k';
hold on
x = 1;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'ko');
h.MarkerFaceColor = [1 1 1];

x = ones(1,sum(~rewMiceInd))*2;
y = msThresh(~rewMiceInd);
h = plot(x,y,'ko');
h.MarkerFaceColor = 'k';
hold on
x = 2;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'ko');
h.MarkerFaceColor = [1 1 1];

figXAxis([],'',[0 3],1:2,{'no rew','rew'})
figYAxis([],'Threshold',[10 24]);
figAxForm([])
title('Visual Trials')

for im = 1:nmice
    subLow(im) = msMatchedHR(im).sub(1);
    subHigh(im) = msMatchedHR(im).sub(2);
end

subplot 222
x = ones(1,sum(rewMiceInd));
y = subLow(rewMiceInd);
h = plot(x,y,'ko');
h.MarkerFaceColor = 'k';
hold on
x = 1;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'ko');
h.MarkerFaceColor = [1 1 1];

x = ones(1,sum(~rewMiceInd))*2;
y = subLow(~rewMiceInd);
h = plot(x,y,'ko');
h.MarkerFaceColor = 'k';
hold on
x = 2;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'ko');
h.MarkerFaceColor = [1 1 1];

x = ones(1,sum(rewMiceInd))*3;
y = subHigh(rewMiceInd);
h = plot(x,y,'bo');
h.MarkerFaceColor = 'b';
hold on
x = 3;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'bo');
h.MarkerFaceColor = [1 1 1];

x = ones(1,sum(~rewMiceInd))*4;
y = subHigh(~rewMiceInd);
h = plot(x,y,'bo');
h.MarkerFaceColor = 'b';
hold on
x = 4;
ybar = mean(y);
yerr = ste(y,2);
h = errorbar(x,ybar,yerr,'bo');
h.MarkerFaceColor = [1 1 1];

figXAxis([],'<<Low,High>>',[0 5],1:4,{'no rew','rew','no rew','rew'})
figYAxis([],'Valid HR - Invalid HR',[-0.2 0.4]);
figAxForm([])
title('Visual Trials')
print('Z:\Analysis\_temp figs\180122\threshold&HR','-dpdf','-fillpage')
%% psychometric across mice 
hit = successIx_all;
miss = missedIx_all;
invHit = invHitIx_all;
invMiss = invMissIx_all;

visTrials = visTargets_all > 0;
audTrials = audTargets_all > 0;
invVisTrials = invVisTargets_all > 0;
invAudTrials = invAudTargets_all > 0;

visBinEdges = exp(linspace(log(min(visTargets)-1),log(max(visTargets)),nBins+1));
audBinEdges = exp(linspace(log(min(audTargets)-(0.5*min(audTargets))),...
    log(max(audTargets)),nBins+1));
[~,~,visBinInd] = histcounts(visTargets_all,visBinEdges);
[~,~,audBinInd] = histcounts(audTargets_all,audBinEdges);
[~,~,invVisBinInd] = histcounts(invVisTargets_all,visBinEdges);
[~,~,invAudBinInd] = histcounts(invAudTargets_all,audBinEdges);

nVisHits = nan(1,nBins);
nVisMisses = nan(1,nBins);
nAudHits = nan(1,nBins);
nAudMisses = nan(1,nBins);
visTargetsBinned = nan(1,nBins);
audTargetsBinned = nan(1,nBins);
visTargetsSte = nan(1,nBins);
audTargetsSte = nan(1,nBins);
nInvVisHits = nan(1,nBins);
nInvVisMisses = nan(1,nBins);
nInvAudHits = nan(1,nBins);
nInvAudMisses = nan(1,nBins);
invVisTargetsBinned = nan(1,nBins);
invAudTargetsBinned = nan(1,nBins);
invVisTargetsSte = nan(1,nBins);
invAudTargetsSte = nan(1,nBins);
for ibin = 1:nBins
    ind = (hit | miss) & visBinInd == ibin;
    if sum(ind) > minTrN_all
    nVisHits(ibin) = sum(hit & visBinInd == ibin);
    nVisMisses(ibin) = sum(miss & visBinInd == ibin);
    visTargetsBinned(ibin) = mean(visTargets_all(ind));
    visTargetsSte(ibin) = ste(visTargets_all(ind),2);
    end
    ind = (hit | miss) & audBinInd == ibin;
    if sum(ind) > minTrN_all
    nAudHits(ibin) = sum(hit & audBinInd == ibin);
    nAudMisses(ibin) = sum(miss & audBinInd == ibin);
    audTargetsBinned(ibin) = mean(audTargets_all(ind));
    audTargetsSte(ibin) = ste(audTargets_all(ind),2);
    end
    ind = (invHit | invMiss) & invVisBinInd == ibin;
    if sum(ind) > minTrN_all
        nInvVisHits(ibin) = sum(invVisBinInd == ibin & invHit);
        nInvVisMisses(ibin) = sum(invVisBinInd == ibin & invMiss);
        invVisTargetsBinned(ibin) = mean(invVisTargets_all(ind));
        invVisTargetsSte(ibin) = ste(invVisTargets_all(ind),2);
    end
    ind = (invHit | invMiss) & invAudBinInd == ibin;
    if sum(ind) > minTrN_all  
        nInvAudHits(ibin) = sum(invAudBinInd == ibin & invHit);
        nInvAudMisses(ibin) = sum(invAudBinInd == ibin & invMiss);
        invAudTargetsBinned(ibin) = mean(invAudTargets_all(ind));
        invAudTargetsSte(ibin) = ste(invAudTargets_all(ind),2);
    end
end

visInd = ~isnan(nVisHits);
invVisInd = ~isnan(nInvVisHits); 
audInd = ~isnan(nAudHits);
invAudInd = ~isnan(nInvAudHits); 

[visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
[audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
[invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
    nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
[invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
    nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));

figure
subplot 121
x = visTargetsBinned(visInd);
xErr = visTargetsSte(visInd);
y = visHR.*100;
yErrL = y - (visHR95ci(:,1)'.*100);
yErrU = (visHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
x = invVisTargetsBinned(invVisInd);
xErr = invVisTargetsSte(invVisInd);
y = invVisHR.*100;
yErrL = y - (invVisHR95ci(:,1)'.*100);
yErrU = (invVisHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
ax = gca;
ax.XScale = 'log';
figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Visual Trials')

subplot 122
x = audTargetsBinned(audInd);
xErr = audTargetsSte(audInd);
y = audHR.*100;
yErrL = y - (audHR95ci(:,1)'.*100);
yErrU = (audHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
x = invAudTargetsBinned(invAudInd);
xErr = invAudTargetsSte(invAudInd);
y = invAudHR.*100;
yErrL = y - (invAudHR95ci(:,1)'.*100);
yErrU = (invAudHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
ax = gca;
ax.XScale = 'log';
figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Auditory Trials')

% print('\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\_temp figs\180118\HRxLevelAllTrials','-dpdf','-fillpage');

% print('\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\_temp figs\180118\HRxStimNAllTrials','-dpdf','-fillpage');

%% psychometric across mice, no short trials;
invLongTrialsInd = invCycles_all > 2;
visTargets = visTargets_all(invLongTrialsInd);
invVisTargets = invVisTargets_all(invLongTrialsInd);
[~,~,visBinInd] = histcounts(visTargets,visBinEdges);
[~,~,invVisBinInd] = histcounts(invVisTargets,visBinEdges);
hitLong = hit(invLongTrialsInd);
missLong = miss(invLongTrialsInd);
invHitLong = invHit(invLongTrialsInd);
invMissLong = invMiss(invLongTrialsInd);

audTargets = audTargets_all(invLongTrialsInd);
invAudTargets = invAudTargets_all(invLongTrialsInd);
[~,~,audBinInd] = histcounts(audTargets,audBinEdges);
[~,~,invAudBinInd] = histcounts(invAudTargets,audBinEdges);

nVisHits = nan(1,nBins);
nVisMisses = nan(1,nBins);
nInvVisHits = nan(1,nBins);
nInvVisMisses = nan(1,nBins);
visTargetsBinned = nan(1,nBins);
visTargetsSte = nan(1,nBins);
invVisTargetsBinned = nan(1,nBins);
invVisTargetsSte = nan(1,nBins);

nAudHits = nan(1,nBins);
nAudMisses = nan(1,nBins);
nInvAudHits = nan(1,nBins);
nInvAudMisses = nan(1,nBins);
audTargetsBinned = nan(1,nBins);
audTargetsSte = nan(1,nBins);
invAudTargetsBinned = nan(1,nBins);
invAudTargetsSte = nan(1,nBins);

for ibin = 1:nBins
    ind = (hitLong | missLong) & visBinInd == ibin;
    if sum(ind) > minTrN_all
        nVisHits(ibin) = sum(hitLong & visBinInd == ibin);
        nVisMisses(ibin) = sum(missLong & visBinInd == ibin);
        visTargetsBinned(ibin) = mean(visTargets(ind));
        visTargetsSte(ibin) = ste(visTargets(ind),2);
    end
    ind = (invHitLong | invMissLong) & invVisBinInd == ibin;
    if sum(ind) > minTrN_all
        nInvVisHits(ibin) = sum(invHitLong & invVisBinInd == ibin);
        nInvVisMisses(ibin) = sum(invMissLong & invVisBinInd == ibin);
        invVisTargetsBinned(ibin) = mean(invVisTargets(ind));
        invVisTargetsSte(ibin) = ste(invVisTargets(ind),2);
    end
    
    ind = (hitLong | missLong) & audBinInd == ibin;
    if sum(ind) > minTrN_all
        nAudHits(ibin) = sum(hitLong & audBinInd == ibin);
        nAudMisses(ibin) = sum(missLong & audBinInd == ibin);
        audTargetsBinned(ibin) = mean(audTargets(ind));
        audTargetsSte(ibin) = ste(audTargets(ind),2);
    end
    ind = (invHitLong | invMissLong) & invAudBinInd == ibin;
    if sum(ind) > minTrN_all
        nInvAudHits(ibin) = sum(invHitLong & invAudBinInd == ibin);
        nInvAudMisses(ibin) = sum(invMissLong & invAudBinInd == ibin);
        invAudTargetsBinned(ibin) = mean(invAudTargets(ind));
        invAudTargetsSte(ibin) = ste(invAudTargets(ind),2);
    end
end
visInd = ~isnan(nVisHits);
invVisInd = ~isnan(nInvVisHits);
audInd = ~isnan(nAudHits);
invAudInd = ~isnan(nInvAudHits);

[visHR,visHR95ci] = binofit(nVisHits(visInd),nVisHits(visInd)+nVisMisses(visInd));
[invVisHR,invVisHR95ci] = binofit(nInvVisHits(invVisInd),...
    nInvVisHits(invVisInd)+nInvVisMisses(invVisInd));
[audHR,audHR95ci] = binofit(nAudHits(audInd),nAudHits(audInd)+nAudMisses(audInd));
[invAudHR,invAudHR95ci] = binofit(nInvAudHits(invAudInd),...
    nInvAudHits(invAudInd)+nInvAudMisses(invAudInd));

figure
suptitle('Excluding short (<= 2 cycles) invalid trials')
subplot 121
x = visTargetsBinned(visInd);
xErr = visTargetsSte(visInd);
y = visHR.*100;
yErrL = y - (visHR95ci(:,1)'.*100);
yErrU = (visHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
x = invVisTargetsBinned(invVisInd);
xErr = invVisTargetsSte(invVisInd);
y = invVisHR.*100;
yErrL = y - (invVisHR95ci(:,1)'.*100);
yErrU = (invVisHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
ax = gca;
ax.XScale = 'log';
figXAxis([],'Orientation Change (deg)',visLevels_lim,visLevels_label,visLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Visual Trials')

subplot 122
x = audTargetsBinned(audInd);
xErr = audTargetsSte(audInd);
y = audHR.*100;
yErrL = y - (audHR95ci(:,1)'.*100);
yErrU = (audHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
x = invAudTargetsBinned(invAudInd);
xErr = invAudTargetsSte(invAudInd);
y = invAudHR.*100;
yErrL = y - (invAudHR95ci(:,1)'.*100);
yErrU = (invAudHR95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,xErr,xErr,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
ax = gca;
ax.XScale = 'log';
figXAxis([],'Tone Volume (?)',audLevels_lim,audLevels_label,audLevels_label)
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label)
figAxForm
title('Auditory Trials')

% print('\\CRASH.dhe.duke.edu\data\home\ashley\Analysis\_temp figs\180118\HRxLevelExcludingShortInv','-dpdf','-fillpage');

%% across mice time bin
cycles = unique(nCycles_all);
nCycles = max(cycles);

visCycHR = nan(1,nCycles);
nVisCycHits = nan(1,nCycles);
nVisCycMisses = nan(1,nCycles);
audCycHR = nan(1,nCycles);
nAudCycHits  = nan(1,nCycles);
nAudCycMisses = nan(1,nCycles);


invVisCycHR = nan(1,nCycles);
nInvVisCycHits = nan(1,nCycles);
nInvVisCycMisses = nan(1,nCycles);
invAudCycHR = nan(1,nCycles);
nInvAudCycHits  = nan(1,nCycles);
nInvAudCycMisses = nan(1,nCycles);

for icyc = 1:nCycles
    nVisCycHits(icyc) = sum(visTrials & nCycles_all == icyc & hit);
    nVisCycMisses(icyc) = sum(visTrials & nCycles_all == icyc & miss);
    ind = visTrials & nCycles_all == icyc & (hit | miss);
    if sum(ind) > minTrN_all
        visCycHR(icyc) = sum(visTrials & nCycles_all == icyc & hit) ./ ...
            sum(visTrials & nCycles_all == icyc & (hit | miss));
    end
    nAudCycHits(icyc) = sum(audTrials & nCycles_all == icyc & hit);
    nAudCycMisses(icyc) = sum(audTrials & nCycles_all == icyc & miss);
    ind = audTrials & nCycles_all == icyc & (hit | miss);
    if sum(ind) > minTrN_all
        audCycHR(icyc) = sum(audTrials & nCycles_all == icyc & hit) ./ ...
            sum(audTrials & nCycles_all == icyc & (hit | miss));
    end

    nInvVisCycHits(icyc) = sum(invVisTrials & invCycles_all == icyc & invHit);
    nInvVisCycMisses(icyc) = sum(invVisTrials & invCycles_all == icyc & invMiss);
    ind = invVisTrials & invCycles_all == icyc & (invHit | invMiss);
    if sum(ind) > minTrN_all
        invVisCycHR(icyc) = sum(invVisTrials & invCycles_all == icyc & invHit) ./ ...
            sum(invVisTrials & invCycles_all == icyc & (invHit | invMiss));
    end
    nInvAudCycHits(icyc) = sum(invAudTrials & invCycles_all == icyc & invHit);
    nInvAudCycMisses(icyc) = sum(invAudTrials & invCycles_all == icyc & invMiss);
    ind = invAudTrials & invCycles_all == icyc & (invHit | invMiss);
    if sum(ind) > minTrN_all
        invAudCycHR(icyc) = sum(invAudTrials & invCycles_all == icyc & invHit) ./ ...
            sum(invAudTrials & invCycles_all == icyc & (invHit | invMiss));
    end
end
[~,visCyc95ci] = binofit(nVisCycHits,nVisCycHits+nVisCycMisses);
[~,audCyc95ci] = binofit(nAudCycHits,nAudCycHits+nAudCycMisses);
[~,invVisCyc95ci] = binofit(nInvVisCycHits,nInvVisCycHits+nInvVisCycMisses);
[~,invAudCyc95ci] = binofit(nInvAudCycHits,nInvAudCycHits+nInvAudCycMisses);


figure;
suptitle(sprintf('Minimum trial number = %s, all levels included',num2str(minTrN_all)))
subplot 121
x = 1:nCycles;
y = visCycHR.*100;
yErrL = y - (visCyc95ci(:,1)'.*100);
yErrU = (visCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
y = invVisCycHR.*100;
yErrL = y - (invVisCyc95ci(:,1)'.*100);
yErrU = (invVisCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],x,x);
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label);
figAxForm
title('Visual Trials')

subplot 122
x = 1:nCycles;
y = audCycHR.*100;
yErrL = y - (audCyc95ci(:,1)'.*100);
yErrU = (audCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
y = invAudCycHR.*100;
yErrL = y - (invAudCyc95ci(:,1)'.*100);
yErrU = (invAudCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],x,x);
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label);
figAxForm
title('Auditory Trials')


%% across mice broad time-bin, level-matched within bin
cycBinEdges = [1 3 4 5 6 7 11];
nBins = length(cycBinEdges);
[~,~,cycBinInd] = histcounts(nCycles_all,cycBinEdges);
[~,~,invCycBinInd] = histcounts(invCycles_all,cycBinEdges);

visCycHR = nan(1,nBins);
nVisCycHits = nan(1,nBins);
nVisCycMisses = nan(1,nBins);
audCycHR = nan(1,nBins);
nAudCycHits  = nan(1,nBins);
nAudCycMisses = nan(1,nBins);


invVisCycHR = nan(1,nBins);
nInvVisCycHits = nan(1,nBins);
nInvVisCycMisses = nan(1,nBins);
invAudCycHR = nan(1,nBins);
nInvAudCycHits  = nan(1,nBins);
nInvAudCycMisses = nan(1,nBins);

for ibin = 1:nBins
    
    nInvVisCycHits(ibin) = sum(invVisTrials & invCycBinInd == ibin & invHit);
    nInvVisCycMisses(ibin) = sum(invVisTrials & invCycBinInd == ibin & invMiss);
    invInd = invVisTrials & invCycBinInd == ibin & (invHit | invMiss);
    if sum(invInd) > minTrN_all
        invVisCycHR(ibin) = sum(invVisTrials & invCycBinInd == ibin & invHit) ./ ...
            sum(invVisTrials & invCycBinInd == ibin & (invHit | invMiss));
    end
    ind = find(visTrials & cycBinInd == ibin & (hit | miss));
    matches = cell2mat(getMatchedValidTrialIndex(visTargets_all(ind),invVisTargets_all(invInd)));
    hitTemp = hit(ind);
    missTemp = miss(ind);
    nVisCycHits(ibin) = sum(hitTemp(matches));
    nVisCycMisses(ibin) = sum(missTemp(matches));
    if length(ind) > minTrN_all
        visCycHR(ibin) = nVisCycHits(ibin) ./ (nVisCycHits(ibin)+nVisCycMisses(ibin));
    end
    
    nInvAudCycHits(ibin) = sum(invAudTrials & invCycBinInd == ibin & invHit);
    nInvAudCycMisses(ibin) = sum(invAudTrials & invCycBinInd == ibin & invMiss);
    invInd = invAudTrials & invCycBinInd == ibin & (invHit | invMiss);
    if sum(invInd) > minTrN_all
        invAudCycHR(ibin) = sum(invAudTrials & invCycBinInd == ibin & invHit) ./ ...
            sum(invAudTrials & invCycBinInd == ibin & (invHit | invMiss));
    end
    ind = find(audTrials & cycBinInd == ibin & (hit | miss));
    matches = cell2mat(getMatchedValidTrialIndex(audTargets_all(ind),invAudTargets_all(invInd)));   
    hitTemp = hit(ind);
    missTemp = miss(ind);
    nAudCycHits(ibin) = sum(hitTemp(matches));
    nAudCycMisses(ibin) = sum(missTemp(matches));
    if sum(ind) > minTrN_all
        audCycHR(ibin) = nAudCycHits(ibin) ./ (nAudCycHits(ibin)+nAudCycMisses(ibin));
    end
end
[~,visCyc95ci] = binofit(nVisCycHits,nVisCycHits+nVisCycMisses);
[~,audCyc95ci] = binofit(nAudCycHits,nAudCycHits+nAudCycMisses);
[~,invVisCyc95ci] = binofit(nInvVisCycHits,nInvVisCycHits+nInvVisCycMisses);
[~,invAudCyc95ci] = binofit(nInvAudCycHits,nInvAudCycHits+nInvAudCycMisses);


figure;
suptitle(sprintf('Minimum trial number = %s, all levels included',num2str(minTrN_all)))
subplot 121
x = cycBinEdges;
y = visCycHR.*100;
yErrL = y - (visCyc95ci(:,1)'.*100);
yErrU = (visCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
y = invVisCycHR.*100;
yErrL = y - (invVisCyc95ci(:,1)'.*100);
yErrU = (invVisCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],x,x);
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label);
figAxForm
title('Visual Trials')

subplot 122
x = cycBinEdges;
y = audCycHR.*100;
yErrL = y - (audCyc95ci(:,1)'.*100);
yErrU = (audCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{val};
h.MarkerFaceColor = [1 1 1];
hold on
y = invAudCycHR.*100;
yErrL = y - (invAudCyc95ci(:,1)'.*100);
yErrU = (invAudCyc95ci(:,2)'.*100) - y;
h = errorbar(x,y,yErrL,yErrU,'o');
h.Color = cueColors{inv};
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Stimulus Number',[0 nCycles+1],x,x);
figYAxis([],'Hit Rate (%)',HR_lim,HR_label,HR_label);
figAxForm
title('Auditory Trials')