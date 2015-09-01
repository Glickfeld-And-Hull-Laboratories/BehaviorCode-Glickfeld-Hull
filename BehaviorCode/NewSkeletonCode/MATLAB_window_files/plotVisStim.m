function input = plotVisStim(data_struct, input)

% essential consts
figNum = 4;
name = 'VisStim';
consts = visstim_constants;
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
nTrial = length(input.tNStimAccepted);
tTrialN = input.trialSinceReset;
% if ~input.doMovingDots
  az_mat = [input.gratingAzimuthDeg:input.gratingAzimuthStepDeg:input.gratingAzimuthDeg+(input.gratingAzimuthStepDeg*(input.gratingAzimuthStepN-1))];
  el_mat = [input.gratingElevationDeg:input.gratingElevationStepDeg:input.gratingElevationDeg+(input.gratingElevationStepDeg*(input.gratingElevationStepN-1))];
  dir_mat = [input.gratingDirectionDeg:input.gratingDirectionStepDeg:input.gratingDirectionDeg+(input.gratingDirectionStepDeg*(input.gratingDirectionStepN-1))];
  diam_mat = [input.gratingDiameterDeg:input.gratingDiameterStepDeg:input.gratingDiameterDeg+(input.gratingDiameterStepDeg*(input.gratingDiameterStepN-1))];
  sf_mat = [input.gratingSpatialFreqCPD.*(input.gratingSpatialFreqStepLog.^(double(input.gratingSpatialFreqStepDir).*((1:input.gratingSpatialFreqStepN)-1)))];
  tf_mat = [double(input.gratingTemporalFreqCPS).*double(input.gratingTemporalFreqStepLog).^(double(input.gratingTemporalFreqStepDir).*((1:double(input.gratingTemporalFreqStepN))-1))];
  con_mat = [double(input.gratingContrast).*(double(input.gratingContrastStepLog).^(double(input.gratingContrastStepDir).*(0:double(input.gratingContrastStepN)-1)))];

  if input.doRetStim
  azs = az_mat;
  els = el_mat;
  else
  azs = az_mat(1);
  els = el_mat(1);
  end

  if input.doDirStim
  dirs = dir_mat;
  else
  dirs = dir_mat(1);
  end

  
    if input.doSizeStim
  diams = diam_mat;
  else
  diams = diam_mat(1);
  end

  if input.doSFStim
  SFs = sf_mat;
  else
  SFs = sf_mat(1);
  end

  if input.doTFStim
  TFs = tf_mat;
  else
  TFs = tf_mat(1);
  end

  if input.doConStim
  cons = con_mat;
  else
  cons = con_mat(1);
  end

  if input.doMatrix
  ncond = (length(azs)*length(els))*length(dirs)*length(diams)*length(SFs)*length(TFs)*length(cons);
  else  
  ncond = max([(length(azs)*length(els)) (length(dirs)) length(diams) length(SFs) length(TFs) length(cons)]);
  end

% not sure why these variables aren't where they should be
% else
%   az_mat = [input.tDotAzimuthDeg(1)];
%   el_mat = [input.tDotElevationDeg(1)];
%   dir_mat = [input.dotDirectionDeg:input.dotDirectionStepDeg:input.dotDirectionDeg+(input.dotDirectionStepDeg*(input.dotDirectionStepN-1))];
%   diam_mat = [input.dotFieldSizeDeg];
%   speed_mat = [input.dotSpeedDPS.*(input.dotSpeedStepLog.^(double(input.dotSpeedqStepDir).*((1:input.dotSpeedStepN)-1)))];
%   con_mat = [double(input.dotContrast).*(double(input.dotContrastStepLog).^(double(input.dotContrastStepDir).*(0:double(input.dotContrastStepN)-1)))];
%   coh_mat = [double(input.dotCoherence).*(double(input.dotCoherenceStepLog).^(double(input.dotCoherenceStepDir).*(0:double(input.dotCoherenceStepN)-1)))];
%   size_mat = [input.dotSizeDeg:input.dotSizeStepDeg:input.dotSizeDeg+(input.dotSizeStepDeg*(input.dotSizeStepN-1))];
%   dens_mat = [input.dotDensityDeg:input.dotDensityStepDeg:input.dotDensityDeg+(input.dotDensityStepDeg*(input.dotDensityStepN-1))];

%   azs = az_mat(1);
%   els = el_mat(1);
%   diams = diam_mat(1);

%   if input.doDotsDirectionStim
%   dirs = dir_mat;
%   else
%   dirs = dir_mat(1);
%   end

%   if input.doDotsContrastStim
%   cons = con_mat;
%   else
%   cons = con_mat(1);
%   end

%   if input.doDotsSizeStim
%   sizes = size_mat;
%   else
%   sizes = size_mat(1);
%   end

%   if input.doDotsSpeedStim
%   speeds = speed_mat;
%   elseif input.doRunFeedback
%   speeds = 'feedback';
%   else
%   speeds = speeds_mat(1);
%   end

%   if input.doDotsDensityStim
%   dens = dens_mat;
%   else
%   dens = dens_mat(1);
%   end

%   if input.doDotsCoherenceStim
%   cohs = coh_mat;
%   else
%   cohs = coh_mat(1);
%   end

% end



wheelIx = celleqel2mat_padded(input.tStimWheelCounter);
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

        text(0.00, 1.0, {'Subject Number'}, 'FontSize', 12);
	text(0.70, 1.0, ...
             { sprintf('%2d', input.subjectNum)}, 'FontSize', 12);

        %if ~input.doMovingDots
        tStr = sprintf( ['Drifting gratings \n', ...
                 'Frames On/Off: %3.0f, %3.0f \n' ,...
        				 'Position (Az,El): %s , %s \n',...
        				 'Diameter (deg): %s \n', ...
        				 'Direction (deg): %s \n', ...
        				 'Spatial freq (cpd) = %s \n', ...
        				 'Temporal freq (cps) = %s \n', ...
        				 'Contrast = %s \n'], ...
                        input.nScansOn, ...
                        input.nScansOff, ...
                        mat2str(azs), mat2str(els), ...
                        mat2str(diams), ...
                        mat2str(dirs), ...
                        mat2str(SFs), ...
                        mat2str(TFs), ...
                        mat2str(cons));
      % else
      %   tStr = sprintf( ['Moving dots  \n', ...
      %            'Frames On/Off: %3.0f, %3.0f \n' ,...
      %            'Position (Az,El): %s , %s \n',...
      %            'Field Diameter (deg): %s \n', ...
      %            'Dot Size (deg): %s \n', ...
      %            'Dot density: %s \n', ...
      %            'Direction (deg): %s \n', ...
      %            'Speed (dps) = %s \n', ...
      %            'Contrast = %s \n',...
      %            'Coherence = %s \n'], ...
      %                   input.nScansOn, ...
      %                   input.nScansOff, ...
      %                   mat2str(azs), mat2str(els), ...
      %                   mat2str(diams), ...
      %                   mat2str(sizes), ...
      %                   mat2str(dens), ...
      %                   mat2str(dirs), ...
      %                   mat2str(speeds), ...
      %                   mat2str(cons),...
      %                   mat2str(cohs));
      % end

        text(0, 0.8, tStr, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);
         
        set(gcf, 'Visible', 'on');

end			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific plots for specific tasks in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 - smoothed perf curve
axH = subplot(spSz{:},2);
hold on;
plot(smooth(double(wheelIx), ceil(nTrial/10), smoothType));
lH = plot(smooth(double(wheelIx), nTrial, smoothType));
set(lH, 'Color', 'r', ...
        'LineWidth', 3);
lH2 = plot(smooth(double(wheelIx), 100, smoothType));
set(lH2, 'Color', 'k', ...
        'LineWidth', 2);
        
title('Running rate by trial')
ylabel('Rate');
xlabel('Trials')        



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
          changedStrOut{end+1} = sprintf('Trial %3d: %s: -> %s', ...
                                         tDesN, tVarName, ...
                                         tChangedToStr);
      end
    end
    if trPer80Changed && tDesN > 1
      tStr = mat2str(input.trPer80History{tDesN});
      changedStrOut{end+1} = sprintf('Tr %3d - trPer80LevelN -> %s', ...
                                     tDesN, tStr);
    end
    if block2TrPer80Changed && tDesN > 1
      tStr = mat2str(input.block2TrPer80History{tDesN});
      changedStrOut{end+1} = sprintf('Tr %3d - Block2TrPer80LevelN -> %s', ...
                                     tDesN, tStr);
    end
  end
  input.changedStr = changedStrOut;
  text(0, 1.01, 'Changed Variables List:', 'FontSize', 12);
  text(0, .95,input.changedStr, ...
       'HorizontalAlignment', 'left', ...
       'VerticalAlignment', 'top');
end


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



