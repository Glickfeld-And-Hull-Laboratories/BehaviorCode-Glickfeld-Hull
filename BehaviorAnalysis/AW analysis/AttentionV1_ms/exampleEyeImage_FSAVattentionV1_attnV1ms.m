clear all
close all
useRandSeed = false;
ds = 'FSAV_attentionV1';
eval(ds)
rc = behavConstsAV;
imgParams_FSAV

titleStr = ds(6:end);
fnin = fullfile(rc.caOutputDir, ds,[titleStr '_eye_']);
load([fnin 'eyeStruct'])

fnout = fullfile(rc.ashley,'Manuscripts','Attention V1','Matlab Figs');
if useRandSeed
    load(fullfile(fnout,'eyeStats'))
    rng(eyeStats.randGeneratorSeed);
end
%%
calib = 1/26.6; %mm per pix
exFrame = 58860;
exExptN = 11;
exMouse = expt(exExptN).mouse;
exDate = expt(exExptN).date;
exRun = expt(exExptN).runs(1,:);
radRange = expt(exExptN).eyeradrange;
%%

load(fullfile(rc.ashleyData,exMouse,'two-photon imaging',exDate,exRun,[exRun '_000_000_eye.mat']))
load(fullfile(rc.ashleyAnalysis,exMouse,'two-photon imaging',exDate,exRun,'eyeTC.mat'))

data = squeeze(data);
[center,radii,metric] = imfindcircles(squeeze(data(:,:,exFrame)),radRange);

% xc = size(data,2)/2;       % image center
% yc = size(data,1)/2;
% W=40;
% offset = W;
% data2 = data(yc-W:yc+W,xc-W:xc+W,:);

exImage = data(:,:,exFrame);
figure;imagesc(exImage); colormap gray
viscircles(center,radii,'EdgeColor','w');
print(fullfile(fnout,'eyeImageCircle'),'-dpdf')

[ypix,xpix] = size(exImage);

sb_xpix_2mm = round(2/calib);
sb_ypix = 10;
sb_xstart = 30;
sb_ystart = ypix - 30;

sbImage = zeros(ypix,xpix);
sbImage(sb_ystart:(sb_ystart+sb_ypix-1),sb_xstart:(sb_xstart+sb_xpix_1mm-1)) = 1;
figure;imagesc(sbImage); colormap gray
print(fullfile(fnout,'eyeImage_scalebar_2mm'),'-dpdf')

%%
writetiff(exImage,fullfile(fnout,'eyeImage'))
writetiff(sbImage,fullfile(fnout,'eyeImage_scalebar_2mm'))