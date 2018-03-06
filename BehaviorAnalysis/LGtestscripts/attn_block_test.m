% three prob
% mouse = 'i550';
% %date = strvcat('161004', '161005', '161006', '161007', '161009', '161013', '161014', '161016', '161017', '161018');
% %time = strvcat('1453', '1430', '1509', '1330', '1135', '1517', '1305', '1304', '1218', '1311');
% dates = strvcat('161018', '161020', '161021', '161023', '161024', '161025', '161026', '161027', '161028', '161110', '161111', '161113', '161115');
% time = strvcat('1311', '1426', '1308','1248', '1140', '1316', '1153', '1430', '1411', '1301', '1211', '1219', '1421');
% trials = {['all'], ['all'],['all'], ['all'], ['all'], 1:325, ['all'], ['all'], 1:250, ['all'], ['all'], ['all'], ['all']};

% mouse = 'i538';
% dates = strvcat('161127','161128','161130', '161201', '161202', '161204', '161205', '161207', '161208', '161209', '161211');
% time = strvcat('1227', '1311', '1316', '1409', '1406', '1620', '1340', '1214', '1250', '1357','1329');
% trials = {1:320, 23:250, 1:270, 1:310,1:325, 1:420, 1:280, 1:225, 1:250, 1:325, 1:325};

% mouse = 'i553';
% dates = strvcat('161202','161205', '161206', '161207', '161208', '161209', '161211');
% time = strvcat('1326', '1222', '1357', '1321', '1357', '1249', '1222');
% trials = {1:225, 1:425, 1:425, ['all'], 1:475, 50:400, ['all']};

%two prob
% mouse = 'i553';
% dates = strvcat('161212','161213', '161214', '161215','161216', '161218','161219','161220');
% time = strvcat('1351', '1208', '1416', '1317', '1232', '1335', '1248', '1356');
% trials = {1:500, 1:475, 1:200, 1:350, ['all'],['all'], ['all'], ['all']};

% mouse = 'i538';
% dates = strvcat('161216', '161218','161219');
% time = strvcat('1342', '1439', '1344');
% trials = {1:400, 23:250, 1:280};

%preblock
% mouse = 'i553';
% dates = strvcat('161220','161221','161223','161226','161227','161228','161229','161230','170101','170102','170103');
% time = strvcat('1356','1338','1117','1107','1041','1041','1246','1207','1200','comb','1304');
% trials = {['all'],['all'],['all'],['all'],['all'],['all'],['all'],['all'],['all'],['all'],'1:300'};

mouse = 'i553';
dates = strvcat('170126','170127','170131','170201','170202','170203','170205');
time = strvcat('1301','1237','1321','1230','1232','1255','1239');
trials = {['all'],['all'],'1:420',['all'],'1:400','1:300',['all']};



clear temp
base = '\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data';
for id = 1:size(dates,1)
    ls = dir(fullfile(base, ['data-' mouse '-' dates(id,:) '-*']));
    if length(ls)==1
        load(fullfile(base,ls.name));
    elseif length(ls)>1
        if time(id) == 'comb';
            for i = 1:length(ls)
                load(fullfile(base,ls(i).name))
                temp(i) = input;
            end
            input = concatenateDataBlocks(temp);
        else
            for i = 1:length(ls)
                if time(id) == ls(i).name(18:21)
                    load(fullfile(base,ls(i).name))
                end
            end
        end
    elseif length(ls) == 0
        error(['no files from' dates(id)])
    end
            
    if id>1
        [a b c d] = fieldComp(input,temp);
        ind = size(c,2);
        if ind>0
            for i = 1:ind
                if isfield(input, c{1,i})
                    input = rmfield(input,c{1,i});
                end
            end
        end
    end
    if  ~ischar(trials{id})
        input = trialChopper(input,trials{id});
    end
    temp(id) = input;
end
input = concatenateDataBlocks(temp);

SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'incorrect');
IIx = strcmp(input.trialOutcomeCell, 'ignore');
probs = cell2mat(input.ProbList);
nprob = length(probs);
tprob = celleqel2mat_padded(input.tStimProbAvgLeft);
tcon = chop(celleqel2mat_padded(input.tGratingContrast),2);
tleft = celleqel2mat_padded(input.tLeftTrial);
leftResp = zeros(size(tleft));
leftResp(intersect(find(SIx),find(tleft))) = 1;
leftResp(intersect(find(MIx),find(tleft == 0))) = 1;
actcon = tcon;
actcon(find(tleft)) = -1.*(tcon(find(tleft)));
cons = unique(actcon);
ncon = length(cons);

allprob = zeros(nprob,ncon,2);
rightprob = zeros(nprob,ncon);
rightci = zeros(ncon,2,nprob);
ntrials = zeros(ncon,nprob);
col_mat = strvcat('b','k','g');
figure;
for iprob = 1:nprob
    ind_p = find(tprob == probs(iprob));
    for icon = 1:ncon
        ind_c = find(actcon == cons(icon));
        ind_pc = intersect(find(MIx+SIx),intersect(ind_p,ind_c));
        ntrials(icon,iprob) = sum(SIx(ind_pc),2)+sum(MIx(ind_pc),2);
        nleft = sum(leftResp(ind_pc),2);
        allprob(iprob,icon,:) = [ntrials(icon,iprob) nleft];
    end
    [rightprob(iprob,:), rightci(:,:,iprob)] = binofit(allprob(iprob,:,1)-allprob(iprob,:,2), allprob(iprob,:,1));
    errorbar(cons, rightprob(iprob,:), rightprob(iprob,:)'-rightci(:,1,iprob), rightci(:,2,iprob)-rightprob(iprob,:)',['-o' col_mat(iprob)])
    hold on
end
ylabel('Probability Right Choice')
legend([num2str((100-probs*100)') repmat('% right', [nprob 1])])
title([mouse ': ' dates(1,:) '-' dates(end,:) '- ' num2str(sum(allprob(:,:,1),2)') ' trials' ])
print(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Behavior\Attn\' date '_' mouse 'probRight.pdf'], '-dpdf')

figure; 
abscon = unique(abs(cons));
nconabs = length(abscon);
allprob = zeros(nprob, nconabs,2,2);
corrprob = zeros(nprob,nconabs,2);
corrci = zeros(nprob,2,nconabs,2);
for icon = 1:nconabs
    subplot(1,nconabs,icon)
    ind_c = find(tcon == abscon(icon));
    for iprob = 1:nprob
        ind_p = find(tprob == probs(iprob));
        ind_pc = intersect(find(MIx+SIx),intersect(ind_p,ind_c));
        for iside = 1:2
            if iside == 1
                ind_s = find(actcon<0);
            else
                ind_s = find(actcon>0);
            end
            ind_pcs = intersect(ind_pc,ind_s);
            ntrials = sum(SIx(ind_pcs),2)+sum(MIx(ind_pcs),2);
            ncorr = sum(SIx(ind_pcs),2);
            allprob(iprob,icon,iside,:) = [ntrials ncorr];
            [corrprob(iprob,icon,iside), corrci(iprob,:,icon,iside)] = binofit(allprob(iprob,icon,iside,2), allprob(iprob,icon,iside,1));
        end
        errorbar([-30 30], corrprob(iprob,icon,:), squeeze(corrprob(iprob,icon,:))-squeeze(corrci(iprob,1,icon,:)), squeeze(corrci(iprob,1,icon,:))-squeeze(corrprob(iprob,icon,:)),['-o' col_mat(iprob)])
        hold on
    end
    title([num2str(chop(abscon(icon),2)*100) '% Con'])
    ylim([0 1])
end
subplot(1,nconabs,1)
ylabel('Probability Correct Choice')
suptitle([mouse ': ' dates(1,:) '-' dates(end,:)])
legend([num2str((100-probs*100)') repmat('% right', [nprob 1])],'Location', 'SouthEast')
print(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Behavior\Attn\' date '_' mouse 'CorrByCon.pdf'], '-dpdf')



        
