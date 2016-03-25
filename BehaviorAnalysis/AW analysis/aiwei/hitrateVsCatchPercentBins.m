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
    measured= input.measuredDegs;

    
    val = max(measured);
    edges = [0 val/3 2*val/3 val+1];
    
    deg = max(uniquedegrees);
    edgesDeg = [0 deg/2 deg+1];
    [M, edgesB] = histc(uniquedegrees, edgesDeg)
    
    
    color = ['k','b','g','r','c','m','y'];
    dotcolor =  ['ok','ob','og','or','oc','om','oy'];
    totalN = zeros(1,7)
    
    degrees1 = uniquedegrees(find(edgesB==1))
    degrees2 = uniquedegrees(find(edgesB==2))
    
    Muniquedegrees{imouse} = uniquedegrees;
    Mhitrates{imouse} = hitrates;
    Mcounts{imouse} = counts;
    Mmeasured{imouse} = measured;
    Mdegrees1{imouse}= degrees1;
    Mdegrees2{imouse}= degrees2;
    
    
figure;

       
%         index = find(celleqel2mat_padded(input.tCatchGratingDirectionDeg)==deg);
      
        
     
        index = [];
        for i= 1:length(degrees1)
            deg = degrees1(i)
            A = find(input.uniqueDeg == deg)  
            index = [index A]
        end
        index = sort(index)
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
            semA(ibin) = std(a(:,ind),[],2)./sqrt(length(ind));
            avglevelone(ibin) = mean(HR(ind),2);
            stdlevelone(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
        end 
        
        labels = cellstr(num2str(N(1:3)'))
        
        
       
        errorbarxy(avgA, avglevelone, semA, stdlevelone,{'ok', 'k', 'k'});
        hold on
        text(avgA, avglevelone, labels, 'horizontal','left', 'vertical','bottom')
        hold on

  
     ylabel('Hit Rate %')
    xlabel('Actual % of catch')
  
    xlim([0 val])
    
    hold on
        clear index
        
        index = [];
        for i= 1:length(degrees2)
            deg = degrees2(i)
            B = find(input.uniqueDeg == deg)  
            index = [index B]
        end
        index = sort(index)

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
            semA(ibin) = std(a(:,ind),[],2)./sqrt(length(ind));
            avglevelone(ibin) = mean(HR(ind),2);
            stdlevelone(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
        end 
        
        labels = cellstr(num2str(N(1:3)'))
        
      
       
        errorbarxy(avgA, avglevelone, semA, stdlevelone,{'ob', 'b', 'b'});
        hold on
        text(avgA, avglevelone, labels, 'horizontal','left', 'vertical','bottom')
        hold on

  
     ylabel('Hit Rate %')
    xlabel('Actual % of catch')
    
    xlim([0 val])
    legend('0-45 degrees', '45-90 degrees')
    title(['i' num2str(mouse_name) '- HR vs Catch % '])
end 

    
    Muniquedegrees{imouse} = uniquedegrees;
    Mhitrates{imouse} = hitrates;
    Mcounts{imouse} = counts;
    Mmeasured{imouse} = measured;
    Mdegrees1{imouse}= degrees1;
    Mdegrees2{imouse}= degrees2;
    
    uniquedegrees_all= cat(2, cell2mat(Muniquedegrees));
    hitrates_all = cat(2, cell2mat(Mhitrates));
    counts_all= cat(2, cell2mat(Mcounts));
    measured_all= cat(2, cell2mat(Mmeasured));
    degrees1_all= cat(2, cell2mat(Mdegrees1));
    degrees2_all = cat(2, cell2mat(Mdegrees2));
    
    val = max(measured_all);
    edges = [0 val/3 2*val/3 val+1];
    
    deg = max(uniquedegrees_all);
    edgesDeg = [0 deg/2 deg+1];
    [M, edgesB] = histc(uniquedegrees_all, edgesDeg)
    
    
    color = ['k','b','g','r','c','m','y'];
    dotcolor =  ['ok','ob','og','or','oc','om','oy'];
    totalN = zeros(1,7)
    
    figure;
    
    index = [];
        for i= 1:length(degrees1_all)
            deg = degrees1_all(i)
            A = find(uniquedegrees_all == deg)  
            index = [index A]
        end
    index = sort(index)
    HR= hitrates_all(index); 
    a = measured_all(index);
        
    bins = [1 2 3];
        
    [N,edgesA] = histc(a,edges);
    avgA_all = zeros(1,max(bins,[],2));
    semA_all = zeros(1,max(bins,[],2));
    avglevelone_all= zeros(1,max(bins,[],2));
    stdlevelone_all= zeros(1,max(bins,[],2));
    
        
       
        for ibin = 1:max(edgesA,[],2)
            
            ind = find(edgesA == ibin)
            avgA_all(ibin) = mean(a(:,ind),2);
            semA_all(ibin) = std(a(:,ind),[],2)./sqrt(length(ind));
            avglevelone_all(ibin) = mean(HR(ind),2);
            stdlevelone_all(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
        end 
        
        labels = cellstr(num2str(N(1:3)'))
        
        
       
        errorbarxy(avgA_all, avglevelone_all, semA_all, stdlevelone_all,{'ok', 'k', 'k'});
        hold on
        text(avgA_all, avglevelone_all, labels, 'horizontal','left', 'vertical','bottom')
        hold on

  
     ylabel('Hit Rate %')
    xlabel('Actual % of catch')
  
    xlim([0 val])
    
    hold on
        clear index
        
        index = [];
        for i= 1:length(degrees2_all)
            deg = degrees2_all(i)
            B = find(uniquedegrees_all == deg)  
            index = [index B]
        end
        index = sort(index)

        HR= hitrates_all(index); 
        a = measured_all(index);
        
        bins = [1 2 3];
        
        [N,edgesA] = histc(a,edges);
        avgA_all = zeros(1,max(bins,[],2));
        semA_all = zeros(1,max(bins,[],2));
        avglevelone_all= zeros(1,max(bins,[],2));
        stdlevelone_all= zeros(1,max(bins,[],2));
    
        
        for ibin = 1:max(edgesA,[],2)
            
            ind = find(edgesA == ibin)
            avgA_all(ibin) = mean(a(:,ind),2);
            semA_all(ibin) = std(a(:,ind),[],2)./sqrt(length(ind));
            avglevelone_all(ibin) = mean(HR(ind),2);
            stdlevelone_all(ibin) = std(HR(ind),[],2)/sqrt(length(HR(ind)));
        end 
        
        labels = cellstr(num2str(N(1:3)'))
        
      
       
        errorbarxy(avgA_all, avglevelone_all, semA_all, stdlevelone_all,{'ob', 'b', 'b'});
        hold on
        text(avgA_all, avglevelone_all, labels, 'horizontal','left', 'vertical','bottom')
        hold on
    
    xlim([0 val])
    legend('0-45 degrees', '45-90 degrees')
    title(['All mice- HR vs Catch % '])





