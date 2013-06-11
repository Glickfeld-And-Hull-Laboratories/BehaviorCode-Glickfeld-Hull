function str = formatBlock2ParamsFromIndex(x1d, b2Num)

if b2Num == 1, 
    exT = 'b2One'; 
else
    exT = 'b2Two'; 
end

idStr = x1d.([exT 'RampOrTrain']);
switch idStr
    case 'ramp'
        switch x1d.([exT 'RampDoExpRamp'])
            case 1
                rStr = 'exp';
            case 0
                rStr = 'lin';
            otherwise
                rStr = 'n/a';
        end

        str = sprintf('ramp: %d (%s)+%d, u/d %d', ...
            x1d.([exT 'RampLengthMs']), ...
            rStr, ...
            x1d.([exT 'RampExtraConstantLengthMs']), ...
            x1d.rampTrainUpDownTimeMs);
    
    case 'train'
        str = sprintf('train: %d pulses (%d ms) period %d ms, u/d %d', ...
            x1d.([exT 'TrainNPulses']), ...
            x1d.([exT 'TrainPulseLengthMs']), ...
            x1d.([exT 'TrainPeriodMs']), ...
            x1d.rampTrainUpDownTimeMs);
end
