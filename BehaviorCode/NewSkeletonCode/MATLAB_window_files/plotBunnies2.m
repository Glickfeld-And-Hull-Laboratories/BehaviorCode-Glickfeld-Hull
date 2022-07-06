function input = plotVisStim(data_struct, input)

% essential consts
figNum = 4;
name = 'TwoStim';
%consts = visstim_constants;

tHostname = lower(hostname);
tUsername = getenv('USER');
homeDir = fullfile('/Users',tUsername);
constsCentralDataPath = fullfile(homeDir, 'Documents/MWorks');
constsDataPath = fullfile(homeDir, 'Documents/MWorks/Data');
constsBehavPdfPath = fullfile(homeDir, 'Documents/MWorks/BunniesOutputPdfs');

spSz = {2,2};

smoothType = 'lowess';
%% draw figure
figH = figure(figNum);
set(figH, 'ToolBar', 'none');
clf;
figPos = [700 140 780 770];
set(figH, 'Position', figPos);
orange = [1 0.64 0];

%% set up arrays
if ~isfield(input, 'changedStr')
input.changedStr = {};
end

% some data processing
nTrial = length(input.ntrials);
tTrialN = input.trialSinceReset;

cItiStart_mat = unique(celleqel2mat_padded(input.cItiStart));
cStimOneOn_mat = unique(celleqel2mat_padded(input.cStimOneOn));
cStimOneOff_mat = unique(celleqel2mat_padded(input.cStimOneOff));
cStimTwoOn_mat = unique(celleqel2mat_padded(input.cStimTwoOn));
cStimTwoOff_mat = unique(celleqel2mat_padded(input.cStimTwoOff));
tisiTimeMs_mat = unique(celleqel2mat_padded(input.tisiTimeMs));
tstimOne_mat = unique(celleqel2mat_padded(input.tstimOne));
tstimTwo_mat = unique(celleqel2mat_padded(input.tstimTwo));
tisiTimeFrames_mat = unique(celleqel2mat_padded(input.tisiTimeFrames));
titiTimeFrames_mat = unique(celleqel2mat_padded(input.titiTimeFrames));
tstimOnTimeFrames_mat = unique(celleqel2mat_padded(input.tstimOnTimeFrames));

% if input.doWheelSpeed
% for itrial = 1:nTrial
% wheelIx(itrial) =input.wheelSpeedValues{itrial}(end)-input.wheelSpeedValues{itrial}(1);
% end
% else
% wheelIx = [];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Values

if tTrialN > 1

axH = subplot(spSz{:}, 1);                        % default axes are 0 to 1

set(axH, 'Visible', 'off');
set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

numTrials = nTrial;

%set(gcf, 'Visible', 'off'); % hide figure during text
% drawing - kludge
text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 18);
text(0.70, 1.25, date, 'FontWeight', 'light', 'FontSize', 18);

text(0.00, 1.0, {'Subject Number'}, 'FontSize', 12);
text(0.70, 1.0, ...
     { sprintf('%2d', input.subjectNum)}, 'FontSize', 12);
text(0.00, 0.8, {'Trial Number'}, 'FontSize', 12);
text(0.7, 0.8, ...
     { sprintf('%2d', size(celleqel2mat_padded(input.ntrials),2))}, 'FontSize', 12);

tStr1 = ['Stim One:'];
if 1
tStr1 = sprintf( [tStr1 '\n',...
                  'Isi Time(Ms): %s \n', ...
                  'Isi Time(Frames): %s \n', ...
                  'Iti Time(Frames): %s \n', ...
                  'Stim One: %s \n', ...
                  'Stim Two: %s \n'], ...
                mat2str(tisiTimeMs_mat), ...
                mat2str(tisiTimeFrames_mat), ...
                mat2str(titiTimeFrames_mat),...
                mat2str(tstimOne_mat),...
                mat2str(tstimTwo_mat));
end
text(0, 0.6, tStr1, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 - smoothed perf curve
% display(cell2mat(input.tisiTimeMs))
% for i = 1:length(input.tstimOne)
%     input.tstimOne{1, i} = cast(input.tstimOne{1, i}, 'double');
% end
% for i = 1:length(input.tstimOne)
%     input.tstimTwo{1, i} = cast(input.tstimTwo{1, i}, 'double');
% end


axH = subplot(spSz{:},2);
hold on;
plot(1:1:length(cell2mat(input.tstimOne)), cell2mat(input.tstimOne));
title('Stim One Presented Per Trial');
ylabel('Stim One');
xlabel('Trial');

axH = subplot(spSz{:},3);
hold on;
plot(1:1:length(cell2mat(input.tstimTwo)), cell2mat(input.tstimTwo));
title('Stim Two Presented Per Trial');
ylabel('Stim Two');
xlabel('Trial');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic code below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List changed text
axH = subplot(spSz{:}, 4);
hold on;
set(axH, 'Visible', 'off');
if nTrial > 1
% recompute text list from all changes
desIx = ~cellfun(@isempty, input.constChangedOnTrial);
desIx(1) = false; % skip the first trial
desTrNs = find(desIx);

changedStrOut = {};
for iT = 1:length(desTrNs)
tDesN = desTrNs(iT);
tChangedList = input.constChangedOnTrial{tDesN};
[nRows nCols] = size(tChangedList);

trPer80Changed = false;
block2TrPer80Changed = false;
for iR = 1:nRows
tVarName = tChangedList{iR,1};
tChangedTo = tChangedList{iR,2};
% special case odds changes
if strfind(tVarName, 'trPer80Level')
trPer80Changed = true;  % and iterate; summarize below
elseif strfind(tVarName, 'block2TrPer80Level')
block2TrPer80Changed = true; % and iterate; summarize below
else
if ischar(tChangedTo)
tChangedToStr = sprintf('''%s''', tChangedTo);
else
tChangedToStr = sprintf('%g', tChangedTo);
end
if ~strcmp(tChangedToStr, '''DEFAULT''')
changedStrOut{end+1} = sprintf('Trial %3d: %s: -> %s', ...
                               tDesN, tVarName, ...
                               tChangedToStr);
end
end
end
end
input.changedStr = changedStrOut;
text(0, 1.01, 'Changed Variables List:', 'FontSize', 12);
text(0, .95,input.changedStr, ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'top');
end


%% Add a save button
if ~exist(constsBehavPdfPath, 'file')
mkdir(constsBehavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  constsBehavPdfPath, ...
                  datestr(now, 'yymmdd'), input.subjectNum);
epParams = { figNum, outName, ...
    'FileFormat', 'pdf', ...
    'Size', [12 12], ...
    'PrintUI', false };
bH = uicontrol(figNum, 'Style', 'pushbutton', ...
               'String', sprintf ('Save PDF figure : %s', outName), ...
               'Units', 'pixels', ...
               'Position', [5 5 650 20], ...
               'Callback', { @saveButtonCb, epParams });

%%%%%%%%%%%%%%%%
%set(gcf, 'Visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subfunctions

function saveButtonCb(hObject, eventdata, epParamsIn)
exportfig_print(epParamsIn{:});

function outM = subCell2PadVect(cellV, padChar, nMaxTrials)
if nargin < 2, padChar = NaN; end
if nargin < 3, nMaxTrials = length(cellV); end
isEIx = cellfun(@isempty, cellV);
outM = repmat(padChar, [1 nMaxTrials]);
outM(~isEIx) = cell2mat(cellV(~isEIx));

% %%%%%%%%%%%%%%%%
% 
% function ebH = subManualErrorbar(lH, upperLength, lowerLength)
% 
% if nargin < 3, lowerLength = upperLength; end
% xv = get(lH, 'XData');
% yv = get(lH, 'YData');
% xvals = cat(1, xv, xv, xv);
% yvals = cat(1, yv-lowerLength(:)', yv+upperLength(:)', yv*NaN);
% xvals = reshape(xvals, [1 prod(size(xvals))]);
% yvals = reshape(yvals, [1 prod(size(yvals))]);
% ebH = plot(xvals, yvals);
% set(ebH, 'Color', get(lH, 'Color'));
% return
% 
% %%%%%%%%%%%%%%%%
% 
% function outStr = subPpGrating(input, doPpBlock2Grating)
% if nargin < 2, doPpBlock2Grating = false; end
% 
% % assemble vars- block2 or not
% tV = struct('gratingDiameterDeg', []);
% if doPpBlock2Grating
% tV.stimOneGratingDirectionDeg = input.block2StimOneGratingDirectionDeg;
% tV.stimOneGratingDiameterDeg = input.block2StimOneGratingDiameterDeg;
% tV.stimOneGratingAzimuthDeg = input.block2StimOneGratingAzimuthDeg;
% tV.stimOneGratingElevationDeg = input.block2StimOneGratingElevationDeg;
% tV.stimOneGratingSpatialFreqCPD = input.block2StimOneGratingSpatialFreqCPD;
% tV.stimOneGratingDurationMs = input.block2StimOneGratingOnTimeMs;
% tV.stimOneGratingTemporalFreqCPS = input.block2StimOneGratingTemporalFreqCPS;
% tV.stimTwoGratingDirectionDeg = input.block2StimTwoGratingDirectionDeg;
% tV.stimTwoGratingDiameterDeg = input.block2StimTwoGratingDiameterDeg;
% tV.stimTwoGratingAzimuthDeg = input.block2StimTwoGratingAzimuthDeg;
% tV.stimTwoGratingElevationDeg = input.block2StimTwoGratingElevationDeg;
% tV.stimTwoGratingSpatialFreqCPD = input.block2StimTwoGratingSpatialFreqCPD;
% tV.stimTwoGratingDurationMs = input.block2StimTwoGratingOnTimeMs;
% tV.stimTwoGratingTemporalFreqCPS = input.block2StimTwoGratingTemporalFreqCPS;
% else
% tV.stimOneGratingDirectionDeg = input.stimOneGratingDirectionDeg;
% tV.stimOneGratingDiameterDeg = input.stimOneGratingDiameterDeg;
% tV.stimOneGratingAzimuthDeg = input.stimOneGratingAzimuthDeg;
% tV.stimOneGratingElevationDeg = input.stimOneGratingElevationDeg;
% tV.stimOneGratingSpatialFreqCPD = input.stimOneGratingSpatialFreqCPD;
% tV.stimOneGratingDurationMs = input.stimOneGratingOnTimeMs;
% tV.stimOneGratingTemporalFreqCPS = input.stimOneGratingTemporalFreqCPS;
% tV.stimTwoGratingDirectionDeg = input.stimTwoGratingDirectionDeg;
% tV.stimTwoGratingDiameterDeg = input.stimTwoGratingDiameterDeg;
% tV.stimTwoGratingAzimuthDeg = input.stimTwoGratingAzimuthDeg;
% tV.stimTwoGratingElevationDeg = input.stimTwoGratingElevationDeg;
% tV.stimTwoGratingSpatialFreqCPD = input.stimTwoGratingSpatialFreqCPD;
% tV.stimTwoGratingDurationMs = input.stimTwoGratingOnTimeMs;
% tV.stimTwoGratingTemporalFreqCPS = input.stimTwoGratingTemporalFreqCPS;
% end
% 
% outStr = sprintf('Stim 1: %gms, %gdeg, %gcpd; Stim 2: %gms, %gdeg, %gcpd' , ...
%                  chop(tV.stimOneGratingDurationMs,2), ...
%                  tV.stimOneGratingDirectionDeg, ...
%                  tV.stimOneGratingSpatialFreqCPD,...
%                  chop(tV.stimTwoGratingDurationMs,2), ...
%                  tV.stimTwoGratingDirectionDeg, ...
%                  tV.stimTwoGratingSpatialFreqCPD);
% 
% %%%%%%%%%%%%%%%%
% 
% function outStr = subPpLaser(input)
% 
% if input.laserRampDoExpRamp == 0,
% rampShapeStr = 'lin';
% else
% rampShapeStr = 'exp';
% end
% 
% if input.laserDoLinearRamp
% outStr = [ sprintf('ramp %d (%s)+ const %d ms', ...
%                    input.laserRampLengthMs, ...
%                    rampShapeStr, ...
%                    input.laserRampExtraConstantLengthMs) ];
% elseif input.laserDoPulseTrain
% outStr = [ sprintf('train, pulse %d, period %d, dur %d', ...
%                    input.laserPulseLengthMs, input.laserPulsePeriodMs, ...
%                    input.laserTrainLengthMs) ];
% end
% 
% 
% %%%%%%%%%%%%%%%%
% 
% function outStr = subPpTrialLaser(input, block2Num)
% 
% tV(1).trialLaserPowerMw = [];
% if block2Num == 1
% prefix = 'trialLaser';
% elseif block2Num == 2
% prefix = 'block2TrialLaser';
% else
% error('bug 1a');
% end
% tV.trialLaserPowerMw = input.([prefix 'PowerMw']);
% tV.trialLaserOnTimeMs = input.([prefix 'OnTimeMs']);
% tV.trialLaserOffTimeMs = input.([prefix 'OffTimeMs']);
% 
% outStr = sprintf('%gmW, on %d off %d', ...
%                  chop(tV.trialLaserPowerMw, 2), ...
%                  chop(tV.trialLaserOnTimeMs, 2), ...
%                  chop(tV.trialLaserOffTimeMs, 2));




