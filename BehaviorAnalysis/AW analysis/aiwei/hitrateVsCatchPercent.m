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

% mouse = [613, 614, 616, 626]
% if length(miceAnalyzed) == 1
%     x = find(mouse == miceAnalyzed);
%     imouse = x;
% else
%     imouse = 1:length(miceAnalyzed);
% end   




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
    
    uniquedegrees = unique(input.uniqueDeg);
    hitrates = input.hitrates;
    counts = input.countDegs;
    measured= input.measuredDegs;

    
    val = max(measured);
    edges = [0 val/3 2*val/3 val+1];
    color = ['k','b','g','r','c','m','y'];
    dotcolor =  ['ok','ob','og','or','oc','om','oy'];
    totalN = zeros(1,7)
    
figure;
    for level = 2: length(uniquedegrees)
        deg = uniquedegrees(level);
%         index = find(celleqel2mat_padded(input.tCatchGratingDirectionDeg)==deg);
        index = find(input.uniqueDeg== deg);
        HR= hitrates(index); 
        a = measured(index);
        
        bins = [1 2 3];
        
        [N,edgesA] = histc(a,edges);
        avgA = zeros(1,max(bins,[],2));
        semA = zeros(1,max(bins,[],2));
        avglevelone= zeros(1,max(bins,[],2));
        stdlevelone= zeros(1,max(bins,[],2));
    
        
        
        
       
        for ibin = 1:max(edgesA,[],2)
            
            ind = find(edgesA == ibin)
            avgA(ibin) = mean(a(:,ind),2);
            semA(ibin) = std(a,[],2)./sqrt(length(ind));
            avglevelone(ibin) = mean(HR(ind),2);
            stdlevelone(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
        end 
        
        labels = cellstr(num2str(N(1:3)'))
        
        c = color(level)
       
        errorbarxy(avgA, avglevelone, semA, stdlevelone,{['o' c], c, c});
        hold on
        text(avgA, avglevelone, labels, 'horizontal','left', 'vertical','bottom')
        hold on
    end
  
     ylabel('Hit Rate %')
    xlabel('Actual % of catch')
    
    xlim([0 val])
    
end   






