function input = plotPavlov(data_struct, input)

% essential consts
figNum = 6;
name = 'Pavlov Stim';
consts = visstim_constants;
spSz = {3,2};

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
nTrial = length(input.tGratingDirectionDeg);
tTrialN = input.trialSinceReset;
dirs = cell2mat(input.gratingDirections);
ddirs = cell2mat(input.dGratingDirections);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performance Values

if tTrialN > 1

        axH = subplot(spSz{:}, 1);						% default axes are 0 to 1
        
	set(axH, 'Visible', 'off');
        set(axH, 'OuterPosition', [0.02 0.75, 0.25, 0.2])

        numTrials = nTrial;
        
        set(gcf, 'Visible', 'off'); % hide figure during text
                                    % drawing - kludge
        text(0.00, 1.25, name, 'FontWeight', 'bold', 'FontSize', 18);
        text(0.70, 1.25, date, 'FontWeight', 'light', 'FontSize', 18);

        elMin = round((now - datenum(input.startDateVec)) * 24*60);
        startStr = datestr(input.startDateVec, 'HH:MM');
        text(0.00, 1.05, {'Subject:', 'Start time + elapsed:', 'Reward vol:'});
  text(0.60, 1.05, ...
             { sprintf('%2d', input.subjectNum), ...
               sprintf('%s + %2dm', startStr, elMin), ...
               sprintf('%.1f ms', input.rewardUs./1000)});

        tStr = sprintf( ['Target grating \n', ...
                 'Stim time: %3.0f \n' ,...
                 'Direction (deg): %s \n', ...
        				 'Position (Az,El): %d , %d \n',...
        				 'Diameter (deg): %d \n', ...
        				 'Spatial freq (cpd) = %d \n', ...
        				 'Speed (dps) = %d \n', ...
        				 'Contrast = %d \n',...
                         'Distractor grating \n', ...
                 'Stim time: %3.0f \n' ,...
                 'Direction (deg): %s \n', ...
        				 'Position (Az,El): %d , %d \n',...
        				 'Diameter (deg): %d \n', ...
        				 'Spatial freq (cpd) = %d \n', ...
        				 'Speed (dps) = %d \n', ...
        				 'Contrast = %d \n'], ...
                        input.gratingDurationMs, ...
                        mat2str(dirs), ...
                        input.gratingAzimuthDeg, input.gratingElevationDeg, ...
                        input.gratingDiameterDeg, ...
                        input.gratingSpatialFreqCPD, ...
                        input.gratingSpeedDPS, ...
                        input.gratingContrast,...
                        input.gratingDurationMs, ...
                        mat2str(ddirs), ...
                        input.gratingAzimuthDeg, input.gratingElevationDeg, ...
                        input.gratingDiameterDeg, ...
                        input.gratingSpatialFreqCPD, ...
                        input.gratingSpeedDPS, ...
                        input.gratingContrast);

        text(0, 0.8, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left');
         
        set(gcf, 'Visible', 'on');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%licks per trial
ndir = length(dirs);
tDirs = celleqel2mat_padded(input.tGratingDirectionDeg);
for idir = 1:ndir
    axH = subplot(3,2,2);
    hold on;
    if nTrial > 1 
        ind = find(tDirs == dirs(idir));
        if length(ind>1)
            nLicks = zeros(length(ind),1);
            for i = 1:length(ind)
                ii = ind(i);
                lickTimes = bsxfun(@minus, double(input.lickometerTimesUs{ii}), double(input.stimOnUs{ii}));
                postStim_licks = find(lickTimes>0);
                preStimEnd_licks = find(lickTimes<input.gratingDurationMs*1000);
                nLicks(i) = length(intersect(postStim_licks,preStimEnd_licks));
            end
            plot(ind, smooth(nLicks,5))
            set(axH, 'Visible', 'on');
            xlabel('Trial Number')
            ylabel('Licks per trial')
        end
    end
end          

%lick rates by stimulus Direction
for idir = 1:ndir
    axH = subplot(3,2, idir+2);
    hold on;
    if nTrial > 1 
        ind = find(tDirs == dirs(idir));
      if length(ind>1)
          lickTimes = [];
          trialNumber = [];
          for i = 1:length(ind)
              ii = ind(i);
              lickTimes = [lickTimes bsxfun(@minus, double(input.lickometerTimesUs{ii}), double(input.stimOnUs{ii}))];
              trialNumber = [trialNumber ii.*ones(1,size(double(input.lickometerTimesUs{ii}),2))];
          end
          if length(lickTimes>1)
            g = gray(nTrial/50);
            for i = 1:20
                if nTrial > 50*i
                    ind_i = find(trialNumber < 50*i);
                   h(i+1) = cdfplot(lickTimes(:,ind_i)./1000);
                   set(h(i+1), 'Color', g(i,:));
                   hold on
                end
            end
            h(1) = cdfplot(lickTimes./1000);
            set(h(1), 'Color', 'r');    
            set(axH, 'Visible', 'on');
            xlim([-2000 4000])
            xlabel('Time from stimulus (ms)')
            ylabel('Fraction of licks')
            title([num2str(dirs(idir)) ' deg'])
          end
        end       
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic code below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% List changed text



%% Add a save button
if ~exist(consts.behavPdfPath, 'file')
  mkdir(consts.behavPdfPath)
end
outName = sprintf('%s/%s-behav-i%03d.pdf', ...
                  consts.behavPdfPath, ...
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

%%%%%%%%%%%%%%%%

function ebH = subManualErrorbar(lH, upperLength, lowerLength)

if nargin < 3, lowerLength = upperLength; end
xv = get(lH, 'XData');
yv = get(lH, 'YData');
xvals = cat(1, xv, xv, xv);
yvals = cat(1, yv-lowerLength(:)', yv+upperLength(:)', yv*NaN);
xvals = reshape(xvals, [1 prod(size(xvals))]);
yvals = reshape(yvals, [1 prod(size(yvals))]);
ebH = plot(xvals, yvals);
set(ebH, 'Color', get(lH, 'Color'));
return

%%%%%%%%%%%%%%%%

function outStr = subPpGrating(input, doPpBlock2Grating)
if nargin < 2, doPpBlock2Grating = false; end

% assemble vars- block2 or not
tV = struct('gratingWidthDeg', []);
if doPpBlock2Grating
  tV.gratingWidthDeg = input.block2GratingWidthDeg;
  tV.gratingHeightDeg = input.block2GratingHeightDeg;
  tV.gratingAzimuthDeg = input.block2GratingAzimuthDeg;
  tV.gratingElevationDeg = input.block2GratingElevationDeg;
  tV.gratingSpatialFreqCPD = input.block2GratingSpatialFreqCPD;
  tV.gratingDurationMs = input.block2GratingDurationMs;
  tV.gratingSpeedDPS = input.block2GratingSpeedDPS;
else
  tV.gratingWidthDeg = input.gratingWidthDeg;
  tV.gratingHeightDeg = input.gratingHeightDeg;
  tV.gratingAzimuthDeg = input.gratingAzimuthDeg;
  tV.gratingElevationDeg = input.gratingElevationDeg;
  tV.gratingSpatialFreqCPD = input.gratingSpatialFreqCPD;
  tV.gratingDurationMs = input.gratingDurationMs;
  tV.gratingSpeedDPS = input.gratingSpeedDPS;
end

outStr = sprintf('%gx%gdeg, at (%g,%g), %gcpd, %dms', ...
                 roundn(double(tV.gratingWidthDeg),-2), ...
                 roundn(double(tV.gratingHeightDeg),-2), ...
                 tV.gratingAzimuthDeg, ...
                 tV.gratingElevationDeg, ...
                 tV.gratingSpatialFreqCPD, ...
                 tV.gratingDurationMs);


if tV.gratingSpeedDPS ~= 0 
  outStr = strcat(outStr, sprintf(' %d deg/s', tV.gratingSpeedDPS));
end

%%%%%%%%%%%%%%%%

function outStr = subPpLaser(input)

if input.laserRampDoExpRamp == 0, 
  rampShapeStr = 'lin'; 
else
  rampShapeStr = 'exp';
end

if input.laserDoLinearRamp
  outStr = [ sprintf('ramp %d (%s)+ const %d ms', ...
                         input.laserRampLengthMs, ...
                         rampShapeStr, ...
                         input.laserRampExtraConstantLengthMs) ];
elseif input.laserDoPulseTrain
  outStr = [ sprintf('train, pulse %d, period %d, dur %d', ...
                         input.laserPulseLengthMs, input.laserPulsePeriodMs, ...
                         input.laserTrainLengthMs) ];
end


%%%%%%%%%%%%%%%%

function outStr = subPpTrialLaser(input, block2Num)

tV(1).trialLaserPowerMw = [];
if block2Num == 1
  prefix = 'trialLaser';
elseif block2Num == 2
  prefix = 'block2TrialLaser';
else
  error('bug 1a');
end
tV.trialLaserPowerMw = input.([prefix 'PowerMw']);
tV.trialLaserOnTimeMs = input.([prefix 'OnTimeMs']);
tV.trialLaserOffTimeMs = input.([prefix 'OffTimeMs']);

outStr = sprintf('%gmW, on %d off %d', ...
                 chop(tV.trialLaserPowerMw, 2), ...
                 chop(tV.trialLaserOnTimeMs, 2), ...
                 chop(tV.trialLaserOffTimeMs, 2));



