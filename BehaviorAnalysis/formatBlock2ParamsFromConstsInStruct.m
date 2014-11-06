function str = formatBlock2ParamsFromConstsInStruct(ds, b2Num)

assert(isempty(b2Num) || b2Num == 1 || b2Num == 2, 'invalid b2Num');

% compat: early b2 vision
if ~isfield(ds, 'doBlock2')
    if isfield(ds, 'block2DoGratingAppearance')
        if ds.block2DoGratingAppearance == 1
            ds.doBlock2 = 1;
        elseif ds.block2DoGratingAppearance == 2
            ds.doBlock2 = 0;
        else
            error('bug');
        end
    end
end


if (~isfield(ds, 'doBlock2') ... % backward compat
    || ds.doBlock2 == 0)
    if(b2Num ~= 1)
        error('MHS:invalidArgument', 'cannot specify b2Num==2 -- b2 not set in data file');
    end

    % compat: some files in 2011
    if isfield(ds, 'doVisualNotChR2')
        error('MHS:badFile', 'old file format that lacks correct constant info - get from doc file');
    end
        
    if isfield(ds,'doLaserStim')
        if ds.doLaserStim
            if ds.laserDoLinearRamp && ~ds.laserDoPulseTrain
                rampOrTrainStr = 'ramp';
                % backward compat
                doExp = 0;
                extraLenMs = 0;
                if isfield(ds, 'laserRampDoExpRamp') 
                    doExp = ds.laserRampDoExpRamp;
                end
                if isfield(ds, 'laserRampExtraConstantLengthMs')
                    extraLenMs = ds.laserRampExtraConstantLengthMs;
                end

                str = subFormatLaserInfo('ramp', ds.laserRampLengthMs, doExp, extraLenMs, ...
                    [], [], [], ds.laserTransitionRampUpDownMs);
                return
            elseif ds.laserDoPulseTrain
                rampOrTrainStr = 'train';
                str = subFormatLaserInfo('train', [], [], [], ...
                    ds.laserTrainLengthMs, ds.laserPulseLengthMs, ds.laserPulsePeriodMs, ...
                    ds.laserTransitionRampUpDownMs);
                return
            else
                error('neither ramp or train is true');
            end
        end
    elseif isfield(ds,'doVisualStim')
        if ds.doVisualStim
            str = subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
                ds.gratingSpatialFreqCPD, ds.gratingDurationMs);
            return
        end
    elseif isfield(ds,'tQuadrature')
        str = subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
            ds.gratingSpatialFreqCPD, ds.reactionTimeMs);
        return
    end
    
elseif isfield(ds, 'doBlock2')
    if ds.doBlock2 == 1
        fList = { 'block2DoGratingAppearance', 'block2DoRampLength', 'block2DoRampVTrain', ...
            'block2DoTrialLaser' };
        for iF = 1:length(fList)
            if ~isfield(ds, fList{iF}), ds.(fList{iF}) = 0; end
        end

    %     if ((ds.block2DoGratingAppearance+ds.block2DoRampLength+ds.block2DoRampVTrain ...
    %         +ds.block2DoTrialLaser) ~= 1)
    %         error('MHS:badFile', 'more than one b2Do param?');
    %     end
        if isfield(ds,'block2StimOnMs')
           ds.block2GratingDurationMs = ds.block2StimOnMs; 
        else
            ds.block2GratingDurationMs = ds.gratingDurationMs;
        end

        if ds.block2DoGratingAppearance == 1 & ds.block2DoTrialLaser==0
            if b2Num == 1
                str = subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
                    ds.gratingSpatialFreqCPD, ds.gratingDurationMs);
            elseif b2Num == 2
                str = subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
                    ds.gratingSpatialFreqCPD, ds.block2GratingDurationMs);
                return
            end
        elseif ds.block2DoRampLength
            error('implement this now');
        elseif ds.block2DoRampVTrain
            if b2Num == 1
                str = subFormatLaserInfo('ramp', ds.laserRampLengthMs, ds.laserRampDoExpRamp, ...
                    ds.laserRampExtraConstantLengthMs, ...
                    [], [], [], ds.laserTransitionRampUpDownMs);
                return
            elseif b2Num == 2
                str = subFormatLaserInfo('train', [], [], [], ...
                    ds.laserTrainLengthMs, ds.laserPulseLengthMs, ds.laserPulsePeriodMs, ...
                    ds.laserTransitionRampUpDownMs);
                return
            end
        elseif ds.block2DoTrialLaser==1 & ds.block2DoGratingAppearance == 0

            tV = struct('trialLaserPowerMw', []);
            fieldList = { 'PowerMw', 'OnTimeMs', 'OffTimeMs' };
            blockPrefix = { 'trialLaser', 'block2TrialLaser'};
            for iB = 1:2
                for iF = 1:length(fieldList)
                    tV(iB).(fieldList{iF}) = ds.([blockPrefix{iB} fieldList{iF}]);
                end
            end

            str = subFormatTrialLaserInfo(tV(b2Num));
        elseif ds.block2DoTrialLaser & ds.block2DoGratingAppearance 

            tV = struct('trialLaserPowerMw', []);
            fieldList = { 'PowerMw', 'OnTimeMs', 'OffTimeMs' };
            blockPrefix = { 'trialLaser', 'block2TrialLaser'};
            for iB = 1:2
                for iF = 1:length(fieldList)
                    tV(iB).(fieldList{iF}) = ds.([blockPrefix{iB} fieldList{iF}]);
                end
            end

            if b2Num == 1
                str = [subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
                    ds.gratingSpatialFreqCPD, ds.gratingDurationMs) [ ] subFormatTrialLaserInfo(tV(b2Num))];
            elseif b2Num == 2
                str = [subFormatVisualInfo(ds.gratingHeightDeg, ds.gratingWidthDeg, ...
                    ds.gratingSpatialFreqCPD, ds.block2GratingDurationMs) [ ] subFormatTrialLaserInfo(tV(b2Num))];
                return
            end

        else
            error('bug');
        end
    end
end

%%% sub functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = subFormatVisualInfo(tHeight, tWidth, tCpd, tDuration)
str = sprintf('grating: %g cpd, %dx%d deg, duration %d ms', ...
    chop(tCpd,2), tWidth, tHeight, tDuration);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = subFormatTrialLaserInfo(tV)
str = sprintf('trialLaser: %g mW', ...
    chop(tV.PowerMw, 3));
    
if tV.OnTimeMs > 0 || tV.OffTimeMs > 0
    str = strcat(str, sprintf(', %d ms on, %d ms off', ...
        tV.OnTimeMs, tV.OffTimeMs));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = subFormatLaserInfo(doRampOrTrain, rampLenMs, rampDoExp, ...
    rampExtraConstLenMs, trainLenMs, ...
    trainPulseLenMs, trainPulsePeriodMs, ...
    transitionUpDownMs);

switch doRampOrTrain
    case 'ramp'
        switch rampDoExp
            case 1
                rStr = 'exp';
            case 0
                rStr = 'lin';
            otherwise
                error('hardcode ramp type if we hit this');
        end
        
        str = sprintf('ramp: %d (%s)+%d ms, u/d %d ms', ...
            rampLenMs, rStr, rampExtraConstLenMs, ...
            transitionUpDownMs);
        return
    case 'train'
        nPulses = floor(trainLenMs / (trainPulseLenMs+(transitionUpDownMs*2)));
        if floor(trainLenMs / (trainPulseLenMs)) ~= nPulses
            error('nPulses different w/ and w/out transition time!!!!!');
        end
        str = sprintf('train: %d pulse(s) (%d ms) period %d ms, u/d %d', ...
            nPulses, trainPulseLenMs, trainPulsePeriodMs, transitionUpDownMs);
        return
end




