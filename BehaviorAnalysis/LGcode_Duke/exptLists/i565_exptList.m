mouse = 'i565';

% add dates {'170801-1322';'170802-1309';'170803-1323';'170804-1319';'170806-1341';
% '170807-1332';'170808-1336';'170809-1315';'170810-1348';'170811-1248'}

%Error using  & 
%NaN's cannot be converted to logicals.

%Error in plot2AFCProbRight (line 107)
%                binofit(sum(pow_mat==pows(ipow) & ratio_mat==ratios(irat) & IiX==0
%                & b2Ix==ib-1 & tRight),totR(ib,irat,ipow));

dates = [170801; 170802; 170803; 170804; 170806; 170807; 170808; 170809; 170810; 170811; ];
trials = {[]; [1 175]; []; []; []; []; []; []; []; []   };

pow_use = [];
mat_use{4} = 2;
mat_use{8} = 2;