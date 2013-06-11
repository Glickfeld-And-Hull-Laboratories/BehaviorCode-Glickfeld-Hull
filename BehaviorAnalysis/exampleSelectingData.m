%% setup
rc = behavConstsHADC8;

% note: to override values in the above file, in your own directory create 
% behavConstsHADC8.m containing:
% 
%%%%%%%%
% function rc = behavConstsHADC8
% % load defaults
% rc = runFromPathOnly('behavConstsHADC8');
% rc.indexFilename = 'whateverYouWant.xls'
% % ... etc.
% return
%%%%%%%

%%%%%%%%%%%
% get data
xd = frm_xls2frm(rc.indexFilename);
fitXd = frm_xls2frm(rc.fitOutputFilename);

% get pulse data
pulseIx = (xd.b2OneRampLengthMs == 0 ...
    & xd.b2OneRampExtraConstantLengthMs == 100 ...
    & xd.b2TwoTrainNPulses > 1);
desIx = pulseIx & ~(xd.UseForPlots == 0 | xd.GoodDataToKeep == 0);
desIx = desIx & xd.b2TwoTrainPeriodMs/2 == xd.b2TwoTrainPulseLengthMs;
desIx = desIx & ismember(xd.Subject,subjNum);

% get data for this
desXdNs = find(desIx);
desXd = frm_extractrow(xd, desXdNs);
% find fit data
desFitNs = frm_findrownum(fitXd, ...
    { 'Subject', xd.Subject(desXdNs), ...
    'DateStr', xd.DateStr(desXdNs), ...
    'UniqueId', xd.UniqueId(desXdNs)});
desFitXd = frm_extractrow(fitXd, desFitNs);

%%%%%%%%%%%%
% draw
figH = figure(9); clf


xsA1 = desXd.b2TwoTrainNPulses;
ys1 = desFitXd.Threshold2;

hold on;
lH = plot(xsA1,  ys1, 'k.-');

set(gca, 'YScale', 'log');
ticklabel_exp2float('Y');
box on;

%% export for keynote
clipclip(gcf, 7*[1 0.75], 'pdf');

