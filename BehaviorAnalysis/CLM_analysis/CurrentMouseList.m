% Current working mice, 2 shifts of 4
mouse = [];
%%
% mouse(1).name = 'i484';
% mouse(1).genotype = 'SOM';
% mouse(1).injection = 'GCaMP7F & tdTomato';
% mouse(1).task = 'Plaid';
% mouse(1).rawtrials = 0;
%%
mouse(1).name = 'i788';
mouse(1).genotype = 'PV';
mouse(1).injection = 'none';
mouse(1).task = 'Ori';
mouse(1).rawtrials = 0;
%%
mouse(2).name = 'i789';
mouse(2).genotype = 'PV';
mouse(2).injection = 'none';
mouse(2).task = 'Ori';
mouse(2).rawtrials = 0;
%%
mouse(3).name = 'i790';
mouse(3).genotype = 'PV';
mouse(3).injection = 'none';
mouse(3).task = 'Ori';
mouse(3).rawtrials = 0;
%%
mouse(4).name = 'i794';
mouse(4).genotype = 'PV';
mouse(4).injection = 'none';
mouse(4).task = 'OriDiscrim';
mouse(4).rawtrials = 0;
%%
mouse(5).name = 'i795';
mouse(5).genotype = 'PV';
mouse(5).injection = 'none';
mouse(5).task = 'OriDiscrim';
mouse(5).rawtrials = 0;
%%
mouse(6).name = 'i792';
mouse(6).genotype = 'PV';
mouse(6).injection = 'none';
mouse(6).task = 'OriDiscrim';
mouse(6).rawtrials = 0;
%%
mouse(7).name = 'i791';
mouse(7).genotype = 'PV';
mouse(7).injection = 'none';
mouse(7).task = 'OriDiscrim';
mouse(7).rawtrials = 0;
%%
mouse(8).name = 'i793';
mouse(8).genotype = 'PV';
mouse(8).injection = 'none';
mouse(8).task = 'OriDiscrim';
mouse(8).rawtrials = 0;
%%
allmice = mouse;
save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\BehaviorSummaries\'...
         'CurrMouseList.mat'],'allmice');