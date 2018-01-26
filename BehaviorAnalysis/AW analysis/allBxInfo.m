clear all
close all

ms2analyze = {'699','730'};
nmice = length(ms2analyze);

rc = behavConstsAV;
ms = struct;
for im = 1:nmice

mouseName = ms2analyze{im};

mworksFilesDir = dir(fullfile(rc.behavData,sprintf('data-i%s-*',mouseName)));
nexp = size(mworksFilesDir,1);
mworksDate = cell(1,nexp);
mworksTime = cell(1,nexp);
catchDay = false(1,nexp);
pctEarly = nan(1,nexp);
pctCorr = nan(1,nexp);
for iexp = 1:nexp    
    
    d = mworksFilesDir(iexp).name(end-14:end-9);
    t = mworksFilesDir(iexp).name(end-7:end-4);
    try
        mworks = loadMworksFile(mouseName,d,t);
    catch
        continue
    end
            
    
    if isfield(mworks,'doShortCatchTrial')
        if mworks.doShortCatchTrial == 1
            mworksDate{iexp} = d;
            mworksTime{iexp} = t;
            catchDay(iexp) = true;
            nt = length(mworks.trialOutcomeCell);
            hit = strcmp(mworks.trialOutcomeCell,'success');
            miss = strcmp(mworks.trialOutcomeCell,'ignore');
            fa = strcmp(mworks.trialOutcomeCell,'failure');
            
            pctEarly(iexp) = sum(fa)./nt;
            pctCorr(iexp) = sum(hit)./nt;            
        end
    end
end

ms(im).mouseName = mouseName;
ms(im).dates = mworksDate;
ms(im).times = mworksTime;
ms(im).catchDay = catchDay;
ms(im).pctEarly = pctEarly;
ms(im).pctCorr = pctCorr;
end