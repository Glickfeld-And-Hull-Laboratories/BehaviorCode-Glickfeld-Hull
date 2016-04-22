 earlyLT50 = 1;
lapseLT10 = 1;
amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
miceAnalyzed = unique(xd.Subject);
av = behavParamsAV;
mice = unique(xd.Subject);
nMice = length(mice);
fn = fullfile(rc.fitOutputSummary, [date '_CatchSummary.mat']);
% fn = fullfile(rc.fitOutputSummary, ['11-Oct-2015_i613_i614_CatchSummary.mat']);
load(fn)

close all

maxEarlyRate =0.5;
maxLapseRate = 0.7;

% mouse = [613, 614, 616, 626]
% if length(miceAnalyzed) == 1
%     x = find(mouse == miceAnalyzed);
%     imouse = x;
% else
%     imouse = 1:length(miceAnalyzed);
% end   

fullgraph = figure;   
for imouse = 1:nMice;
    mouse_name = mice(imouse);
    if earlyLT50
        early_ind = find(mouse(imouse).early_mat<maxEarlyRate);
    else
        early_ind = 1:size(early_mat,1);
    end
    if lapseLT10
        lapse_ind = intersect(find(mouse(imouse).HR_ori_mat>maxLapseRate), find(mouse(imouse).HR_amp_mat>maxLapseRate));
    else
        lapse_ind = 1:size(HR_ori_mat,1);
    end
    
    use_ind = intersect(early_ind,lapse_ind);
    input = mouse(imouse).input(use_ind);
    input = concatenateStructures(input');
    
    uniquedegrees = unique(input.uniqueDeg);
    hitrates = input.hitrates;
    counts = input.countDegs;
    measured= input.measuredcatchvalue;
    missedIx = strcmp(input.trialOutcomeCell, 'ignore');
    successIx = strcmp(input.trialOutcomeCell, 'success');
    earliesIx = strcmp(input.trialOutcomeCell, 'failure');
    FAIx = strcmp(input.catchTrialOutcomeCell, 'FA');
    CRIx = strcmp(input.catchTrialOutcomeCell, 'CR');
    earliesIx(FAIx) = 0;
    earliesIx(CRIx) = 0;
    
    indexes = input.indexes;
    
    catchDirectionDeg = double(celleqel2mat_padded(input.tCatchGratingDirectionDeg));
    catchAmplitude = chop(double(celleqel2mat_padded(input.tSoundCatchAmplitude, NaN, 'double')),2);
    catchOris = unique(catchDirectionDeg);
    catchAmps = unique(catchAmplitude(~isnan(catchAmplitude)));
    targetDirectionDeg = (double(celleqel2mat_padded(input.tGratingDirectionDeg)));
    targetAmplitude = chop(double(celleqel2mat_padded(input.tSoundTargetAmplitude, NaN, 'double')),2);
    targetOris = unique(targetDirectionDeg);
    targetAmps = unique(targetAmplitude(~isnan(targetAmplitude)));
    
    targetReact = celleqel2mat_padded(input.reactTimeMs);
    catchReact = celleqel2mat_padded(input.leverUpTimeMs) - celleqel2mat_padded(input.tCatchTimeMs);
    
    Mmeasured{imouse} = measured;
    
    Indexes{imouse} = indexes;
    TargetDegs = unique(targetDirectionDeg);
    
%%

    
    
% 
%     tarDeg{imouse} = targetDirectionDeg;
%     tarAmp{imouse} = targetAmplitude;
%     cDeg{imouse} = catchDirectionDeg;
%     cAmp{imouse} = catchAmplitude;
%     
%     sIx_all{imouse} = successIx;
%     mIx_all{imouse} = missedIx;
%     eIx_all{imouse} = earliesIx;
%     faIx_all{imouse} = FAIx;
%     crIx_all{imouse} = CRIx;
    
    
  
    
%%    
    
%     measured_all= cat(2, cell2mat(Mmeasured))   
    measured_all = measured;
    val = max(measured);
    edges = [0 val/3 2*val/3 val+1];
    
    bins = [1 2 3];
    
    [N,edgesA] = histc(measured,edges);
            
    avgA = zeros(1,max(bins,[],2));
    semA = zeros(1,max(bins,[],2));
    avglevelone= zeros(1,max(bins,[],2));
    stdlevelone= zeros(1,max(bins,[],2))
    
    binOne = find(measured_all <= edges(2));
    binTwo = find(measured_all <= edges(3) & measured_all > edges(2));
    binThree = find(measured_all <= edges(4) & measured_all > edges(3));
  
    
    TotalFA = [];
    TotalCR= [];
    TotalSuccessOne=[];
    TotalIgnoreOne= [];
    
    TotalTargetReactmean = [];
    TotalTargetReactstd = [];
    TotalCatchReactmean= [];
    TotalCatchReactstd= [];
    
    FAs = [];
    CRs= [];
    SuccessOne = [];
    IgnoreOne= [];
    TRs= [];
    CatchsR = [];

    
    IndexOne = indexes(binOne)
    for x = 1:max(length(IndexOne))
        ind = IndexOne(x)
       
        uniquedegrees = mouse(imouse).input(ind).uniqueDeg;
        catchDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tCatchGratingDirectionDeg);
        
        targetDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tGratingDirectionDeg);
        
        newtargetDirectionDeg = [];
        if length(uniquedegrees) == 0
        end
        
        for i = 1:max(length(uniquedegrees))
            degree = uniquedegrees(i);
            e = find(catchDirectionDeg == degree);
            Catchpercent = length(e) / length(catchDirectionDeg);
            f = find(targetDirectionDeg == degree);
            Targetpercent = length(f) / length(targetDirectionDeg);
            g = randperm(length(f), length(e));
            newtargetDirectionDeg =  [newtargetDirectionDeg f(g) ]  
        end
        
        FAday= sum(mouse(imouse).input(ind).FAIx);
        FAind = find(mouse(imouse).input(ind).FAIx==1);
        
        FAs= [FAs FAday];
        CRday= sum(mouse(imouse).input(ind).CRIx);
        CRs = [CRs CRday];
        Successday = sum(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg));
        SuccessOne = [SuccessOne Successday];
        
        Successind = find(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg)==1);

        Ignoreday = sum(mouse(imouse).input(ind).Ignorex(newtargetDirectionDeg));
        IgnoreOne = [Ignoreday IgnoreOne];
        
        TargetReactday = mouse(imouse).input(ind).reactTimeMs;
        TRs= [TRs TargetReactday(Successind)];
         
        
        CatchReactday = mouse(imouse).input(ind).CatchReactday;
        CatchsR = [CatchsR CatchReactday(FAind)];
    
    
    end
    
    TotalFA = [TotalFA sum(FAs)]
    TotalCR = [TotalCR sum(CRs)]
    TotalSuccessOne = [TotalSuccessOne sum(SuccessOne)]
    TotalIgnoreOne = [TotalIgnoreOne sum(IgnoreOne)]
    
    avgA(1) = mean(measured(binOne));
    semA(1) = std(measured(binOne))./sqrt(length(binOne));
    
    
    TotalTargetReactmean = [TotalTargetReactmean mean(cell2mat_padded(TRs))];
    TotalTargetReactstd = [ TotalTargetReactstd std(cell2mat_padded(TRs)/ sqrt(length(TRs)))];
    TotalCatchReactmean= [ TotalCatchReactmean mean(CatchsR)];
    TotalCatchReactstd= [TotalCatchReactstd std(CatchsR)/ sqrt(length(CatchsR))];    
    
    
    
    FAs = [];
    CRs= [];
    SuccessOne = [];
    IgnoreOne= [];
    TRs= [];
    CatchsR = [];

    clear ind    
    IndexTwo = indexes(binTwo)
    for y = 1:max(length(IndexTwo))
        ind = IndexTwo(y);         
        uniquedegrees = mouse(imouse).input(ind).uniqueDeg;
        catchDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tCatchGratingDirectionDeg);
        
        targetDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tGratingDirectionDeg);
        
        newtargetDirectionDeg = [];
        if uniquedegrees(1) == 0
            uniquedegrees = uniquedegrees(2:length(uniquedegrees))
        end
        
        for i = 1:max(length(uniquedegrees))
            degree = uniquedegrees(i);
            e = find(catchDirectionDeg == degree);
            Catchpercent = length(e) / length(catchDirectionDeg);
            f = find(targetDirectionDeg == degree);
            Targetpercent = length(f) / length(targetDirectionDeg);
            g = randperm(length(f), length(e));
            newtargetDirectionDeg =  [newtargetDirectionDeg f(g) ]  
        end
        
        FAday= sum(mouse(imouse).input(ind).FAIx);
        FAs= [FAs FAday];
        FAind = find(mouse(imouse).input(ind).FAIx==1);
        CRday= sum(mouse(imouse).input(ind).CRIx);
        CRs = [CRs CRday];
        Successday = sum(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg));
        SuccessOne = [SuccessOne Successday];
        Successind = find(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg)==1);

        
        
        Ignoreday = sum(mouse(imouse).input(ind).Ignorex(newtargetDirectionDeg));
        IgnoreOne = [Ignoreday IgnoreOne];
        
        
        TargetReactday = mouse(imouse).input(ind).reactTimeMs;
        TRs= [TRs TargetReactday(Successind)];
        
     
        CatchReactday = mouse(imouse).input(ind).CatchReactday;
        CatchsR = [CatchsR CatchReactday(FAind)];
    end
    TotalFA = [TotalFA sum(FAs)]
    TotalCR = [TotalCR sum(CRs)]
    TotalSuccessOne = [TotalSuccessOne sum(SuccessOne)]
    TotalIgnoreOne = [TotalIgnoreOne sum(IgnoreOne)]
    avgA(2)= mean(measured(binTwo));
    semA(2) = std(measured(binTwo)) ./ sqrt(length(binTwo));
    
    TotalTargetReactmean = [TotalTargetReactmean mean(cell2mat_padded(TRs))];
    TotalTargetReactstd = [TotalTargetReactstd std(cell2mat_padded(TRs)/ sqrt(length(TRs)))];
    TotalCatchReactmean= [ TotalCatchReactmean mean(CatchsR)];
    TotalCatchReactstd= [TotalCatchReactstd std(CatchsR)/ sqrt(length(CatchsR))];    
    
    
    FAs = [];
    CRs= [];
    SuccessOne = [];
    IgnoreOne= [];
    TRs= [];
    CatchsR = [];
    
    
    clear ind    
    IndexThree = indexes(binThree)
    for y = 1:max(length(IndexThree))
        ind = IndexThree(y)
        
        
               
        uniquedegrees = mouse(imouse).input(ind).uniqueDeg;
        catchDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tCatchGratingDirectionDeg);
        
        targetDirectionDeg = celleqel2mat_padded(mouse(imouse).input(ind).tGratingDirectionDeg);
        
        newtargetDirectionDeg = [];
        if uniquedegrees(1) == 0
            uniquedegrees = uniquedegrees(2:length(uniquedegrees))
        end
        
        for i = 1:max(length(uniquedegrees))
            degree = uniquedegrees(i);
            e = find(catchDirectionDeg == degree);
            Catchpercent = length(e) / length(catchDirectionDeg);
            f = find(targetDirectionDeg == degree);
            Targetpercent = length(f) / length(targetDirectionDeg);
            g = randperm(length(f), length(e));
            newtargetDirectionDeg =  [newtargetDirectionDeg f(g) ]  
        end
        
        FAday= sum(mouse(imouse).input(ind).FAIx);
        FAs= [FAs FAday];
        FAind = find(mouse(imouse).input(ind).FAIx==1);
        CRday= sum(mouse(imouse).input(ind).CRIx);
        CRs = [CRs CRday];
        Successday = sum(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg));
        SuccessOne = [SuccessOne Successday];
        Successind = find(mouse(imouse).input(ind).SuccessIx(newtargetDirectionDeg)==1);
        Ignoreday = sum(mouse(imouse).input(ind).Ignorex(newtargetDirectionDeg));
        IgnoreOne = [Ignoreday IgnoreOne];
        
        
                
        TargetReactday = mouse(imouse).input(ind).reactTimeMs;
        TRs= [TRs TargetReactday(Successind)];
        
     
        CatchReactday = mouse(imouse).input(ind).CatchReactday;
        CatchsR = [CatchsR CatchReactday(FAind)];
    end
    TotalFA = [TotalFA sum(FAs)]
    TotalCR = [TotalCR sum(CRs)]
    TotalSuccessOne = [TotalSuccessOne sum(SuccessOne)]
    TotalIgnoreOne = [TotalIgnoreOne sum(IgnoreOne)]
    avgA(3)= mean(measured(binThree));
    semA(3) = std(measured(binThree)) ./ sqrt(length(binThree));
    
    TotalTargetReactmean = [TotalTargetReactmean mean(cell2mat_padded(TRs))];
    TotalTargetReactstd = [TotalTargetReactstd std(cell2mat_padded(TRs)/ sqrt(length(TRs)))];
    TotalCatchReactmean= [ TotalCatchReactmean mean(CatchsR)];
    TotalCatchReactstd= [TotalCatchReactstd std(CatchsR)/ sqrt(length(CatchsR))];    
    
     [HR_ori, ci95_HR_ori] = binofit(TotalSuccessOne, TotalIgnoreOne+ TotalSuccessOne);
    [FR_ori, ci95_FR_ori] = binofit(TotalFA, TotalFA+TotalCR);

figure; 
    errorbarxy(avgA, HR_ori, semA, semA, HR_ori - ci95_HR_ori(:,1)', ci95_HR_ori(:,2)' - HR_ori, {'ok', 'k', 'k'});
    hold on
    errorbarxy(avgA, FR_ori, semA, semA, FR_ori - ci95_FR_ori(:,1)', ci95_FR_ori(:,2)' - FR_ori, {'oc', 'c', 'c'});
  
    title(['i' num2str(mouse_name) '- HR vs Invalid trial % '])
    xlabel('Actual % of Invalid trials')
    ylabel('Hit rate')
    
    ylim([0 1])
    
    

 figure;
  
    errorbarxy(avgA, TotalTargetReactmean, semA, TotalTargetReactstd,{'ok', 'k', 'k'})
    hold on
    errorbarxy(avgA, TotalCatchReactmean, semA, TotalCatchReactstd,{'oc', 'c', 'c'})

        
   title(['i' num2str(mouse_name) '- RT vs Invalid trial % '])
    xlabel('Actual % of Invalid trials')
    ylabel('RT') 
    
    
    
%     for ibin = 1:max(edgesA,[],2)
%         
%         ind = find(edgesA == ibin)
%         avgA(ibin) = mean(a(:,ind),2);
%         semA(ibin) = std(a(:,ind),[],2)./sqrt(length(ind));
%         avglevelone(ibin) = mean(HR(ind),2);
%         stdlevelone(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
%     end
%     
    
    
%     
%       ind_catch = find(n_ori(2,:)>5);
%     HitsV = zeros(1,length(ind_catch));
%     MissesV = zeros(1,length(ind_catch));
%     FAsV = zeros(1,length(ind_catch));
%     CRsV = zeros(1,length(ind_catch));
%     if ~isempty(ind_catch)
%     for ibin = 1:length(ind_catch)
%        tarC = unique(catchDirectionDeg(find(bin_oric==(ind_catch(ibin)))));
%        indT = find(ismember(targetDirectionDeg,tarC));
%        indC = find(ismember(catchDirectionDeg,tarC));
%        HitsV(ibin) = sum(successIx(indT),2);
%        MissesV(ibin) = sum(missedIx(indT),2);
%        FAsV(ibin) = sum(FAIx(indC),2);
%        CRsV(ibin) = sum(CRIx(indC),2); 
%     end
%     [HR_V, HR_ci95_V] = binofit(HitsV, HitsV+MissesV);
%     [FR_V, FR_ci95_V] = binofit(FAsV, FAsV+CRsV);
%     errorbarxy(HR_V, FR_V, HR_V-HR_ci95_V(:,1)', HR_ci95_V(:,2)'-HR_V, FR_V-FR_ci95_V(:,1)', FR_ci95_V(:,2)'-FR_V, {['o' av(ms_ind).col_str], av(ms_ind).col_str, av(ms_ind).col_str})
%     hold on
%     xlim([0 1])
%     ylim([0 1])
%     plot(x,y,'--k')
%     hold on
%     xlabel('Hit rate- valid cue')
%     ylabel('Hit rate- invalid cue')
%     title('Visual trials')
%     end
    
end