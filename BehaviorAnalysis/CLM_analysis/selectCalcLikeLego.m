function [rndSelectivity, rndMaskSelectivity] = selectCalcLikeLego(input)

nTrial = length(input.trialOutcomeCell);
block2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);

block2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
block1Ix = block2Ix==0;

if isfield(input,'isNoGo')
noGoIx = cell2mat_padded(input.isNoGo);
else
noGoIx = zeros(size(block2Ix));
end
zeroConIx = zeros(size(block2Ix));
if input.doZeroConTrials
  totCon = sum(celleqel2mat_padded(input.tGratingContrast)+celleqel2mat_padded(input.dGratingContrast),1);
  zeroConIx(find(totCon==0)) = 1;
end

if isfield(input,'doMask') && input.doMask
    maskIx = celleqel2mat_padded(input.tDoMask);
    type1MaskIx = maskIx;
    if isfield(input,'doType2Mask')
      type2MaskIx = celleqel2mat_padded(input.tDoType2Mask);
      type1MaskIx(find(type2MaskIx)) = 0;
    else
      type2MaskIx = zeros(size(block2Ix));
    end
else
    maskIx = zeros(size(block2Ix));
    type2MaskIx = zeros(size(block2Ix));
    type1MaskIx = zeros(size(block2Ix));
end

correctIx = strcmp(input.trialOutcomeCell, 'success');
incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
missIx = ignoreIx;
nCorr = sum(correctIx);
nInc = sum(incorrectIx);
nIg = sum(ignoreIx);
nTrials = length(input.trialOutcomeCell);
decisionV = celleqel2mat_padded(input.tDecisionTimeMs);
trPer80Str = regexprep(['[' num2str(input.trPer80V, '%3d') ']'], '[ *', '[');
if isfield(input, 'leftTrPer80Level1')
  leftTrPer80Str = ['[' num2str([input.leftTrPer80Level1 input.leftTrPer80Level2 input.leftTrPer80Level3 input.leftTrPer80Level4 input.leftTrPer80Level5 input.leftTrPer80Level6 input.leftTrPer80Level7 input.leftTrPer80Level8]) ']' ];
else
  leftTrPer80Str = [];
end

leftTrialIx = celleqel2mat_padded(input.tLeftTrial);
leftTrialIx(find(noGoIx)) = 0;
leftTrialIx(find(zeroConIx)) = 0;
leftTrNs = find(leftTrialIx);
rightTrialIx = ~leftTrialIx;
rightTrialIx(find(noGoIx)) = 0;
rightTrialIx(find(zeroConIx)) = 0;
rightTrNs = find(rightTrialIx);
leftTrialIx = logical(leftTrialIx);
nLeft = sum(leftTrialIx);
nRight = sum(rightTrialIx);
noGoIx = logical(noGoIx);
nNoGo = sum(noGoIx);

leftOutcomes = input.trialOutcomeCell(leftTrialIx);
leftCorr = strcmp(leftOutcomes, 'success');
nLeftCorr = sum(leftCorr);
leftIgn = strcmp(leftOutcomes, 'ignore');
nLeftIgn = sum(leftIgn);
leftInc = strcmp(leftOutcomes, 'incorrect');
nLeftInc = sum(leftInc);

rightOutcomes = input.trialOutcomeCell(rightTrialIx);
rightCorr = strcmp(rightOutcomes, 'success');
nRightCorr = sum(rightCorr);
rightIgn = strcmp(rightOutcomes, 'ignore');
nRightIgn = sum(rightIgn);
rightInc = strcmp(rightOutcomes, 'incorrect');
nRightInc = sum(rightInc);

correctIx = strcmp(input.trialOutcomeCell, 'success');
incorrectIx = strcmp(input.trialOutcomeCell, 'incorrect');
ignoreIx = strcmp(input.trialOutcomeCell, 'ignore');
missIx = ignoreIx;
nCorr = sum(correctIx);
nInc = sum(incorrectIx);
nIg = sum(ignoreIx);
nTrials = length(input.trialOutcomeCell);



 if input.doContrastDiscrim
    if input.doContrastDetect
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif input.gratingContrastDiffSPO > 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
    elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
      differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
    end
  elseif input.doSizeDiscrim
    if input.gratingDiameterDiffSPO > 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingDiameterDeg) - celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    elseif ~isfield(input, 'dGratingDiameterDiff') & input.gratingDiameterDiffSPO <= 10
      differenceRight = chop(celleqel2mat_padded(input.rightGratingDiameterDeg) - celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    elseif isfield(input, 'dGratingDiameterDiff') & input.gratingDiameterDiffSPO <= 10
        differenceRight =chop(celleqel2mat_padded(input.rightGratingDiameterDeg) ./ celleqel2mat_padded(input.leftGratingDiameterDeg),2);
    end
  elseif input.doOriDiscrim
      differenceRight = chop(celleqel2mat_padded(input.tGratingDirectionStart) - double(input.gratingTargetDirection),2);
      if isfield(input,'doMask')
          differenceRight_mask = chop(celleqel2mat_padded(input.tPlaidDirectionStart) - double(input.gratingTargetDirection),2);
      else
          differenceRight_mask = nan(size(differenceRight));
      end
  else
  if input.gratingContrastDiffSPO > 10
    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
  elseif ~isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
    differenceRight = chop(celleqel2mat_padded(input.rightGratingContrast) - celleqel2mat_padded(input.leftGratingContrast),2);
  elseif isfield(input, 'dGratingContrastDiff') & input.gratingContrastDiffSPO <= 10
    differenceRight =chop(celleqel2mat_padded(input.rightGratingContrast) ./ celleqel2mat_padded(input.leftGratingContrast),2);
  end
end

plotTrsB1 = differenceRight((correctIx|incorrectIx)&~block2Ix&~maskIx);
nLevelsB1 = unique(plotTrsB1);
curr100 = zeros(size(block2Ix));
nT= size(curr100,2);
if nT>100
  curr100(nT-99:nT) = 1;
end
percentContCellB1 = cell(1,length(nLevelsB1));
percentContCellB1_100 = cell(1,length(nLevelsB1));
for kk=1:length(nLevelsB1)
    valB1 = nLevelsB1(kk);
    valIxB1 = differenceRight==valB1;
    totalNTrialsValB1 = length(differenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix&~maskIx));
    totalNTrialsValB1_100 = length(differenceRight(valIxB1&(correctIx|incorrectIx)&~block2Ix&curr100&~maskIx));
    if min(differenceRight) < 0
      if valB1>=0,
          ind = setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]);  
          rightNTrialsValB1 = length(ind);
          ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]));  
          rightNTrialsValB1_100 = length(ind2);
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
          percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
      elseif valB1<0,
          ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]);
          rightNTrialsValB1 = length(ind);
          ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]));
          rightNTrialsValB1_100 = length(ind2);
          percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
          percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
      end
    else
        if valB1>=1,
            ind = setdiff(intersect(find(valIxB1),find(correctIx)), find(block2Ix));
            rightNTrialsValB1 = length(ind);
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
            ind2 = intersect(find(curr100), setdiff(intersect(find(valIxB1),find(correctIx)), [find(block2Ix) find(maskIx)]));
            rightNTrialsValB1_100 = length(ind2);
            percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
        elseif valB1<1,
            ind = setdiff(intersect(find(valIxB1),find(incorrectIx)), find(block2Ix));
            rightNTrialsValB1 = length(ind);
            percentContCellB1{kk} = rightNTrialsValB1/totalNTrialsValB1;
            ind2 = intersect(find(curr100),setdiff(intersect(find(valIxB1),find(incorrectIx)), [find(block2Ix) find(maskIx)]));
            rightNTrialsValB1_100 = length(ind2);
            percentContCellB1_100{kk} = rightNTrialsValB1_100/totalNTrialsValB1_100;
        end
    end
end

if sum(maskIx)>0
    if sum(type1MaskIx)>0
      plotTrsM1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type1MaskIx);
      nLevelsM1 = unique(plotTrsM1);
      percentContCellM1 = cell(1,length(nLevelsM1));
      for kk=1:length(nLevelsM1)
          valM1 = nLevelsM1(kk);
          valIxM1 = differenceRight_mask==valM1;
          totalNTrialsValM1 = length(differenceRight_mask(valIxM1&(correctIx|incorrectIx)&~block2Ix&type1MaskIx));
          if min(differenceRight_mask) < 0
            if valM1>=0,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                rightNTrialsValM1 = length(ind);
                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
            elseif valM1<0,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                rightNTrialsValM1 = length(ind);
                percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
            end
          else
              if valM1>=1,
                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(correctIx))), find(block2Ix));
                  rightNTrialsValM1 = length(ind);
                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
              elseif valM1<1,
                  ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM1),find(incorrectIx))), find(block2Ix));
                  rightNTrialsValM1 = length(ind);
                  percentContCellM1{kk} = rightNTrialsValM1/totalNTrialsValM1;
              end
          end
      end
    end
    if sum(type2MaskIx)>0
      plotTrsT2M1 = differenceRight_mask((correctIx|incorrectIx)&~block2Ix&type2MaskIx);
      nLevelsT2M1 = unique(plotTrsT2M1);
      percentContCellT2M1 = cell(1,length(nLevelsT2M1));
      for kk=1:length(nLevelsT2M1)
          valT2M1 = nLevelsT2M1(kk);
          valIxT2M1 = differenceRight_mask==valT2M1;
          totalNTrialsValT2M1 = length(differenceRight_mask(valIxT2M1&(correctIx|incorrectIx)&~block2Ix&type2MaskIx));
          if min(differenceRight_mask) < 0
            if valT2M1>=0,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                rightNTrialsValT2M1 = length(ind);
                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
            elseif valT2M1<0,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                rightNTrialsValT2M1 = length(ind);
                percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
            end
          else
              if valT2M1>=1,
                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(correctIx))), find(block2Ix));
                  rightNTrialsValT2M1 = length(ind);
                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
              elseif valT2M1<1,
                  ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M1),find(incorrectIx))), find(block2Ix));
                  rightNTrialsValT2M1 = length(ind);
                  percentContCellT2M1{kk} = rightNTrialsValT2M1/totalNTrialsValT2M1;
              end
          end
      end
    end
end

if sum(block2Ix)>0
  plotTrsB2 = differenceRight((correctIx|incorrectIx)&block2Ix);
  nLevelsB2 = unique(plotTrsB2);
  percentContCellB2 = cell(1,length(nLevelsB2));
  for kk=1:length(nLevelsB2)
    valB2 = nLevelsB2(kk);
    valIxB2 = differenceRight==valB2;
    totalNTrialsValB2 = length(differenceRight(valIxB2&(correctIx|incorrectIx)&block2Ix));
    if min(differenceRight) < 0
      if valB2>=0,
          ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
          rightNTrialsValB2 = length(ind);
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      elseif valB2<0,
          ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
          rightNTrialsValB2 = length(ind);
          percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
      end
    else
        if valB2>=1,
            ind = setdiff(intersect(find(valIxB2),find(correctIx)), find(block1Ix));
            rightNTrialsValB2 = length(ind);
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        elseif valB2<1,
            ind = setdiff(intersect(find(valIxB2),find(incorrectIx)), find(block1Ix));
            rightNTrialsValB2 = length(ind);
            percentContCellB2{kk} = rightNTrialsValB2/totalNTrialsValB2;
        end
    end
  end
end

if sum(maskIx&block2Ix)>0
  if sum(type1MaskIx)>0 type1MaskIx
    plotTrsM2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type1MaskIx);
    nLevelsM2 = unique(plotTrsM2);
    percentContCellM2 = cell(1,length(nLevelsM2));
    for kk=1:length(nLevelsM2)
        valM2 = nLevelsM2(kk);
        valIxM2 = differenceRight_mask==valM2;
        totalNTrialsValM2 = length(differenceRight_mask(valIxM2&(correctIx|incorrectIx)&block2Ix&type1MaskIx));
        if min(differenceRight_mask) < 0
          if valM2>=0,
              ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
              rightNTrialsValM2 = length(ind);
              percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
          elseif valM2<0,
              ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
              rightNTrialsValM2 = length(ind);
              percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
          end
        else
            if valM2>=1,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(correctIx))), find(block1Ix));
                rightNTrialsValM2 = length(ind);
                percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
            elseif valM2<1,
                ind = setdiff(intersect(find(type1MaskIx),intersect(find(valIxM2),find(incorrectIx))), find(block1Ix));
                rightNTrialsValM2 = length(ind);
                percentContCellM2{kk} = rightNTrialsValM2/totalNTrialsValM2;
            end
        end
    end
  end
  if sum(type2MaskIx)>0
    plotTrsT2M2 = differenceRight_mask((correctIx|incorrectIx)&block2Ix&type2MaskIx);
    nLevelsT2M2 = unique(plotTrsT2M2);
    percentContCellT2M2 = cell(1,length(nLevelsT2M2));
    for kk=1:length(nLevelsT2M2)
        valT2M2 = nLevelsT2M2(kk);
        valIxT2M2 = differenceRight_mask==valT2M2;
        totalNTrialsValT2M2 = length(differenceRight_mask(valIxT2M2&(correctIx|incorrectIx)&block2Ix&type2MaskIx));
        if min(differenceRight_mask) < 0
          if valT2M2>=0,
              ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
              rightNTrialsValT2M2 = length(ind);
              percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
          elseif valT2M2<0,
              ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
              rightNTrialsValT2M2 = length(ind);
              percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
          end
        else
            if valT2M2>=1,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(correctIx))), find(block1Ix));
                rightNTrialsValT2M2 = length(ind);
                percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
            elseif valT2M2<1,
                ind = setdiff(intersect(find(type2MaskIx),intersect(find(valIxT2M2),find(incorrectIx))), find(block1Ix));
                rightNTrialsValT2M2 = length(ind);
                percentContCellT2M2{kk} = rightNTrialsValT2M2/totalNTrialsValT2M2;
            end
        end
    end
  end
end


didNoGoIx = celleqel2mat_padded(input.didNoGo);
plotTrsNoGo = differenceRight(noGoIx|didNoGoIx);
nLevelsNoGo = unique(plotTrsNoGo);
for kk=1:length(nLevelsNoGo)
  valNoGo = nLevelsNoGo(kk);
  valIx = differenceRight==valNoGo;
  totalNTrialsVal = sum(valIx);
  totalNTrialsValNoGo = sum(valIx & didNoGoIx);
  percentContCellNoGo{kk} = totalNTrialsValNoGo/totalNTrialsVal;
end

maxCell = length(nLevelsB1);
lowMidCell = maxCell/2;
highMidCell = lowMidCell + 1;
A = cell2mat(percentContCellB1(maxCell));
B = cell2mat(percentContCellB1(1));
M = (A + B)/2;
Selectivity = A - B;
rndSelectivity = roundn(Selectivity, -3);
Bias = M - 0.5;
rndBias = roundn(Bias, -2);

if input.doMask && sum(maskIx)>0
    maxCellmask = length(nLevelsM1);
    lowMidCell = maxCellmask/2;
    highMidCell = lowMidCell + 1;
    A = cell2mat(percentContCellM1(maxCellmask));
    B = cell2mat(percentContCellM1(1));
    M = (A + B)/2;
    maskSelectivity = A - B;
    rndMaskSelectivity = roundn(maskSelectivity, -3);
    Bias = M - 0.5;
    rndBias = roundn(Bias, -2);
else
    rndMaskSelectivity = NaN;
end


if mod(highMidCell,1) == 0
  C = cell2mat(percentContCellB1(highMidCell));
  D = cell2mat(percentContCellB1(lowMidCell));
  M2 = (C + D)/2;
  InnerBias = M2 - 0.5;
  rndInnerBias = roundn(InnerBias, -2);
else
  InnerBias = NaN;
  rndInnerBias = NaN;
end
