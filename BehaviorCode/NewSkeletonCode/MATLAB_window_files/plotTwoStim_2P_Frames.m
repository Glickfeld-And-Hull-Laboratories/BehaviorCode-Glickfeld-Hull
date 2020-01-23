function input = plotVisStim(data_struct, input)

% essential consts
figNum = 4;
name = 'TwoStim';
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

  stim1On_mat = unique(celleqel2mat_padded(input.nStimOneFramesOn));
  az1_mat = unique(celleqel2mat_padded(input.tStimOneGratingAzimuthDeg));
  el1_mat = unique(celleqel2mat_padded(input.tStimOneGratingElevationDeg));
  dir1_mat = unique(celleqel2mat_padded(input.tStimOneGratingDirectionDeg));
  diam1_mat = chop(unique(celleqel2mat_padded(input.tStimOneGratingDiameterDeg)),2); 
  sf1_mat = unique(celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD)); 
  tf1_mat = unique(celleqel2mat_padded(input.tStimOneGratingTemporalFreqCPS));
  con1_mat = unique(celleqel2mat_padded(input.tStimOneGratingContrast));
  phase1_mat = unique(celleqel2mat_padded(input.tStimOneGratingPhaseDeg));
  sound1_mat = unique(celleqel2mat_padded(input.tStimOneSoundAmplitude));

  mask1_dir_mat = unique(celleqel2mat_padded(input.tMaskOneGratingDirectionDeg));
  mask1_con_mat = unique(celleqel2mat_padded(input.tMaskOneGratingContrast));
  mask1_phase_mat = unique(celleqel2mat_padded(input.tMaskOneGratingPhaseDeg));

  stim2On_mat = unique(celleqel2mat_padded(input.nStimTwoFramesOn));
  az2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingAzimuthDeg));
  el2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingElevationDeg));
  dir2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingDirectionDeg));
  diam2_mat = chop(unique(celleqel2mat_padded(input.tStimTwoGratingDiameterDeg)),2); 
  sf2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingSpatialFreqCPD)); 
  tf2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingTemporalFreqCPS));
  con2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingContrast));
  phase2_mat = unique(celleqel2mat_padded(input.tStimTwoGratingPhaseDeg));
  sound2_mat = unique(celleqel2mat_padded(input.tStimTwoSoundAmplitude));

  mask2_dir_mat = unique(celleqel2mat_padded(input.tMaskTwoGratingDirectionDeg));
  mask2_con_mat = unique(celleqel2mat_padded(input.tMaskTwoGratingContrast));

  isi_mat = unique(celleqel2mat_padded(input.nFramesISI));

  if input.doWheelSpeed
    for itrial = 1:nTrial
        wheelIx(itrial) =input.wheelSpeedValues{itrial}(end)-input.wheelSpeedValues{itrial}(1);
    end
  else
    wheelIx = [];
  end
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
        text(0.00, 0.8, {'Trial Number'}, 'FontSize', 12);
  text(0.7, 0.8, ...
             { sprintf('%2d', size(celleqel2mat_padded(input.tStimulusNumber),2))}, 'FontSize', 12);

  tStr1 = ['Stim One:'];
        if input.stimOneDoVisualStim
          tStr1 = sprintf( [tStr1 '\n',...
                 'Frames On: %s \n' ,...
        				 'Az/El (deg): %s/%s \n',...
        				 'Diameter (deg): %s \n', ...
        				 'Direction (deg): %s \n', ...
        				 'Spatial freq (cpd): %s \n', ...
        				 'Temporal freq (cps): %s \n', ...
        				 'Contrast (pct): %s \n', ...
                         'Phase (deg): %s \n', ...
                         'Sound Vol: %s'], ...
                        mat2str(stim1On_mat), ...
                        mat2str(az1_mat), mat2str(el1_mat), ...
                        mat2str(diam1_mat), ...
                        mat2str(dir1_mat), ...
                        mat2str(sf1_mat), ...
                        mat2str(tf1_mat), ...
                        mat2str(con1_mat.*100),...
                        mat2str(phase1_mat),...
                        mat2str(sound1_mat));
      elseif input.stimOneDoAuditoryStim
        tStr1 = sprintf( [tStr1 '\n', ...
          'Sound Vol = %s'], ...
            mat2str(sound1_mat));
      else
        tStr1 = sprintf( ['No Stim One \n']);
      end

      tStr2 = ['Stim Two:'];
      if input.stimTwoDoVisualStim
          tStr2 = sprintf( [tStr2 '\n', ... 
                 'Frames On: %s \n',...
                 'Az/El (deg): %s/%s \n',...
                 'Diameter (deg): %s \n', ...
                 'Direction (deg): %s \n', ...
                 'Spatial freq (cpd): %s \n', ...
                 'Temporal freq (cps): %s \n', ...
                 'Contrast (pct): %s \n', ...
                 'Phase (deg): %s \n'], ...
                        mat2str(stim2On_mat), ...
                        mat2str(az2_mat),... 
                        mat2str(el2_mat), ...
                        mat2str(diam2_mat), ...
                        mat2str(dir2_mat), ...
                        mat2str(sf2_mat), ...
                        mat2str(tf2_mat), ...
                        mat2str(con2_mat.*100),...
                        mat2str(phase2_mat),...
                        mat2str(sound2_mat));
      elseif input.stimTwoDoAuditoryStim
        tStr2 = sprintf( ['tStr2 \n', ...
          'Sound Vol = %s'], ...
            mat2str(sound2_mat));
      else
        tStr2 = sprintf( ['No Stim Two \n']);
      end

      tMaskStr1 = ['Mask One:'];
      if input.doMask & find(mask1_con_mat > 0)
          tMaskStr1 = sprintf( [tMaskStr1 '\n',...
                 'Direction (deg): %s \n', ...
                 'Contrast (pct): %s \n', ...
                 'Phase (deg): %s \n'], ...
                        mat2str(mask1_dir_mat), ...
                        mat2str(mask1_con_mat.*100), ...
                        mat2str(mask1_phase_mat));
      else
        tMaskStr1 = sprintf( ['No Mask One \n']);
      end

      tMaskStr2 = ['Mask Two:'];
      if input.doMask & find(mask2_con_mat == 1)
          tMaskStr2 = sprintf( [tMaskStr2 '\n',...
                 'Direction (deg): %s \n', ...
                 'Contrast (pct): %s \n'], ...
                        mat2str(mask2_dir_mat), ...
                        mat2str(mask2_con_mat.*100));
      else
        tMaskStr2 = sprintf( ['No Mask Two \n']);
      end

      tStrISI = sprintf('ISI (frames): %s', mat2str(isi_mat));

        text(0, 0.6, tStr1, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);

        text(0, -.5, tStr2, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);

        text(1, 0.6, tMaskStr1, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);

        text(1, -.5, tMaskStr2, ...
             'VerticalAlignment', 'top', ...
             'HorizontalAlignment', 'left', 'FontSize', 12);

        text(0, -2, tStrISI, ...
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
tV = struct('gratingDiameterDeg', []);
if doPpBlock2Grating
  tV.stimOneGratingDirectionDeg = input.block2StimOneGratingDirectionDeg;
  tV.stimOneGratingDiameterDeg = input.block2StimOneGratingDiameterDeg;
  tV.stimOneGratingAzimuthDeg = input.block2StimOneGratingAzimuthDeg;
  tV.stimOneGratingElevationDeg = input.block2StimOneGratingElevationDeg;
  tV.stimOneGratingSpatialFreqCPD = input.block2StimOneGratingSpatialFreqCPD;
  tV.stimOneGratingDurationMs = input.block2StimOneGratingOnTimeMs;
  tV.stimOneGratingTemporalFreqCPS = input.block2StimOneGratingTemporalFreqCPS;
  tV.stimTwoGratingDirectionDeg = input.block2StimTwoGratingDirectionDeg;
  tV.stimTwoGratingDiameterDeg = input.block2StimTwoGratingDiameterDeg;
  tV.stimTwoGratingAzimuthDeg = input.block2StimTwoGratingAzimuthDeg;
  tV.stimTwoGratingElevationDeg = input.block2StimTwoGratingElevationDeg;
  tV.stimTwoGratingSpatialFreqCPD = input.block2StimTwoGratingSpatialFreqCPD;
  tV.stimTwoGratingDurationMs = input.block2StimTwoGratingOnTimeMs;
  tV.stimTwoGratingTemporalFreqCPS = input.block2StimTwoGratingTemporalFreqCPS;
else
  tV.stimOneGratingDirectionDeg = input.stimOneGratingDirectionDeg;
  tV.stimOneGratingDiameterDeg = input.stimOneGratingDiameterDeg;
  tV.stimOneGratingAzimuthDeg = input.stimOneGratingAzimuthDeg;
  tV.stimOneGratingElevationDeg = input.stimOneGratingElevationDeg;
  tV.stimOneGratingSpatialFreqCPD = input.stimOneGratingSpatialFreqCPD;
  tV.stimOneGratingDurationMs = input.stimOneGratingOnTimeMs;
  tV.stimOneGratingTemporalFreqCPS = input.stimOneGratingTemporalFreqCPS;
  tV.stimTwoGratingDirectionDeg = input.stimTwoGratingDirectionDeg;
  tV.stimTwoGratingDiameterDeg = input.stimTwoGratingDiameterDeg;
  tV.stimTwoGratingAzimuthDeg = input.stimTwoGratingAzimuthDeg;
  tV.stimTwoGratingElevationDeg = input.stimTwoGratingElevationDeg;
  tV.stimTwoGratingSpatialFreqCPD = input.stimTwoGratingSpatialFreqCPD;
  tV.stimTwoGratingDurationMs = input.stimTwoGratingOnTimeMs;
  tV.stimTwoGratingTemporalFreqCPS = input.stimTwoGratingTemporalFreqCPS;
end

outStr = sprintf('Stim 1: %gms, %gdeg, %gcpd; Stim 2: %gms, %gdeg, %gcpd' , ...
                 chop(tV.stimOneGratingDurationMs,2), ...
                 tV.stimOneGratingDirectionDeg, ...
                 tV.stimOneGratingSpatialFreqCPD,...
                 chop(tV.stimTwoGratingDurationMs,2), ...
                 tV.stimTwoGratingDirectionDeg, ...
                 tV.stimTwoGratingSpatialFreqCPD);

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



