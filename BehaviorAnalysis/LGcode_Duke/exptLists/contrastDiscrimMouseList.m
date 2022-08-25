mouse = [];
%% i410 -good - low trials(417)
% mouse(1).name = 'i410';
% mouse(1).genotype = 'SOM';
% mouse(1).opsin = 'ChR2';
% mouse(1).task = 'ContrastDiscrim';
% mouse(1).rawtrials = 417;
%% i414 -good
mouse(1).name = 'i414';
mouse(1).genotype = 'SOM';
mouse(1).opsin = 'ChR2';
mouse(1).task = 'ContrastDiscrim';
mouse(1).rawtrials = 3295;
%% i419 -good
mouse(2).name = 'i419';
mouse(2).genotype = 'PV';
mouse(2).opsin = 'ChR2';
mouse(2).task = 'ContrastDiscrim';
mouse(2).rawtrials = 3301;
%% i441 -good
mouse(3).name = 'i441';
mouse(3).genotype = 'PV';
mouse(3).opsin = 'Arch';
mouse(3).task = 'ContrastDiscrim';
mouse(3).rawtrials = 2748;
%% i442 -good
mouse(4).name = 'i442';
mouse(4).genotype = 'PV';
mouse(4).opsin = 'Arch';
mouse(4).task = 'ContrastDiscrim';
mouse(4).rawtrials = 13906;
%% i443 -good
mouse(5).name = 'i443';
mouse(5).genotype = 'PV';
mouse(5).opsin = 'Arch';
mouse(5).task = 'ContrastDiscrim';
mouse(5).rawtrials = 9359;
%% i459 -good
mouse(6).name = 'i459';
mouse(6).genotype = 'SOM';
mouse(6).opsin = 'Arch';
mouse(6).task = 'ContrastDiscrim';
mouse(6).rawtrials = 3642;
%% i460 -good
mouse(7).name = 'i460';
mouse(7).genotype = 'SOM';
mouse(7).opsin = 'Arch';
mouse(7).task = 'ContrastDiscrim';
mouse(7).rawtrials = 9468;
%% i461 -good
mouse(8).name = 'i461';
mouse(8).genotype = 'SOM';
mouse(8).opsin = 'ChR2';
mouse(8).task = 'ContrastDiscrim';
mouse(8).rawtrials = 2124;
%% i462 -good
mouse(9).name = 'i462';
mouse(9).genotype = 'SOM';
mouse(9).opsin = 'Arch';
mouse(9).task = 'ContrastDiscrim';
mouse(9).rawtrials = 1896;
%% i463 -good
mouse(10).name = 'i463';
mouse(10).genotype = 'SOM';
mouse(10).opsin = 'Arch';
mouse(10).task = 'ContrastDiscrim';
mouse(10).rawtrials = 5032;
%% i464 -dud (1 good day; 229 trials) (DROP)
% mouse(12).name = 'i464';
% mouse(12).genotype = 'SOM';
% mouse(12).opsin = 'ChR2';
% mouse(12).task = 'ContrastDiscrim';
%% i547 - (1 good day; 284 trials)
% mouse(13).name = 'i547';
% mouse(13).genotype = 'PV';
% mouse(13).opsin = 'ChR2';
% mouse(13).task = 'ContrastDiscrim';
%% i548 - good
mouse(11).name = 'i548';
mouse(11).genotype = 'PV';
mouse(11).opsin = 'ChR2';
mouse(11).task = 'ContrastDiscrim';
mouse(11).rawtrials = 1080;
%% i565 - good
mouse(12).name = 'i565';
mouse(12).genotype = 'PV';
mouse(12).opsin = 'ChR2';
mouse(12).task = 'ContrastDiscrim';
mouse(12).rawtrials = 2268;
%% i578 - good
mouse(13).name = 'i578';
mouse(13).genotype = 'PV';
mouse(13).opsin = 'ChR2';
mouse(13).task = 'ContrastDiscrim';
mouse(13).rawtrials = 4300;
%% i581 - good
mouse(14).name = 'i581';
mouse(14).genotype = 'SOM';
mouse(14).opsin = 'ChR2';
mouse(14).task = 'ContrastDiscrim';
mouse(14).rawtrials = 3457;
%% i582 - good
mouse(15).name = 'i582';
mouse(15).genotype = 'SOM';
mouse(15).opsin = 'ChR2';
mouse(15).task = 'ContrastDiscrim';
mouse(15).rawtrials = 3295;
%% i591 -good ; negative control
mouse(16).name = 'i591';
mouse(16).genotype = 'PV';
mouse(16).opsin = 'Arch-';
mouse(16).task = 'ContrastDiscrim';
mouse(16).rawtrials = 3982;
%% i592 - good; negative control
mouse(17).name = 'i592';
mouse(17).genotype = 'PV';
mouse(17).opsin = 'Arch-';
mouse(17).task = 'ContrastDiscrim';
mouse(17).rawtrials = 3000;
%% i593 -good
mouse(18).name = 'i593';
mouse(18).genotype = 'PV';
mouse(18).opsin = 'ChR2';
mouse(18).task = 'ContrastDiscrim';
mouse(18).rawtrials = 2973;
%%
allmice = mouse;
save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\Contrast_Discrim\'...
         'CDMouseList.mat'],'allmice');