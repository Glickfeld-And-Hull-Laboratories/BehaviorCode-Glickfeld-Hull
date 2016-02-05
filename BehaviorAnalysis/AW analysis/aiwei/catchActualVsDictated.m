earlyLT50 = 1;
lapseLT10 = 1;
amp_edges = [0.003 0.008 0.0195 0.047 0.12 0.29 0.73];
ori_edges = [6 12 24 40 60 80 100];

rc = behavConstsAV;
xd = frm_xls2frm(rc.indexFilename, [], rc.indexTextCols);
miceAnalyzed = unique(xd.Subject);
av = behavParamsAV;
fn = fullfile(rc.fitOutputSummary, [date '_i613_i614_CatchSummary.mat']);
% fn = fullfile(rc.fitOutputSummary, ['11-Oct-2015_i613_i614_CatchSummary.mat']);
load(fn)

close all

maxEarlyRate = 1.0;
maxLapseRate = 0.0;

for imouse = 1:length(miceAnalyzed);
    mouse_name = av(imouse).mouse; 
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
    
    
        
  
    catchvalues = input.catchvalue;
    measured= input.measuredcatchvalue;
    whichone = zeros(size(measured));
    
    for idx= 1:length(catchvalues)
        if catchvalues(idx) == 0.025
            whichone(idx) = 1;
        elseif catchvalues(idx)== 0.05 
            whichone(idx) =2;
        elseif catchvalues(idx) == 0.10 
            whichone(idx) = 4;
        end
    end
    
    onex = measured(find(whichone == 1))
    twox = measured(find(whichone == 2))
    fourx= measured(find(whichone == 4))
    
   one = repmat(0.02,1, length(onex))
   two= repmat(0.05, 1, length(twox))
   four = repmat(0.10,1,length(fourx))
   
   onemean = mean(onex)
   onestd = std(onex)/sqrt(length(onex))
   twomean = mean(twox)
   twostd = std(twox)/sqrt(length(twox))
    fourmean = mean(fourx)
   fourstd = std(fourx)/sqrt(length(fourx))
   
   figure;
    scatter(one, onex);
    hold on
    scatter(two, twox);
    hold on 
    scatter(four, fourx);
    hold on
  
    errorbar(0.02, onemean, onestd,onestd,'k*')
    hold on
    errorbar(0.05, twomean, twostd,twostd,'k*')
    hold on
    errorbar(0.10, fourmean, fourstd, fourstd,'k*')
    hold on
    xlim([0 0.10])
    
    x = [0 0.10]
     y = x;
    plot(x,y,'--k')
    hold on
           
            
    
%     index = [1:16]
%     
%     figure;
%     scatter(index, measured);
%     hold on;
%     plot(catchvalues);
%     hold on;
%     ylim([0 0.2])
%     xlim([0 6])
    
end