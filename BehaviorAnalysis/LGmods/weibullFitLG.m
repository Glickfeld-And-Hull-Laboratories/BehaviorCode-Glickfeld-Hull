function outS = weibullFit(intensityV, fractCorr, clampAtZero, use50Thresh, inputWeightsIn);
% fit a weibull; definition below and returned in modelFun.
% extract 50% threshold if use50Thresh == true, else (default) use 63%, p(1)
%
%   outS = weibullFit(intensityV, fractCorr, clampAtZero, use50Thresh); 
%   
%   You will nearly always want to pass {'nTrials', trialVec} as inputWeights -
%   will automatically weight by standard deviation
%
% histed 120528 - new thresh

if nargin<3, clampAtZero = false; end
if nargin<4, use50Thresh = false; end
if nargin <5, inputWeightsIn = ones(size(fractCorr)); end

%% process args
assert(length(intensityV) == length(fractCorr));
assert(all(fractCorr >= 0 & fractCorr <= 1));

% error checking
if length(intensityV) < 3
    error('Must have at least 3 intensity pts');
end
if sum(~isnan(fractCorr)) < 3
    error('Must have at least 3 correct pts');
end

% compute input weights
if isnumeric(inputWeightsIn)
    assert(length(inputWeightsIn) == length(intensityV));
    inputWeights = inputWeightsIn;
elseif iscell(inputWeightsIn)
    assert(length(inputWeightsIn) == 2);
    switch inputWeightsIn{1}
        case 'nTrials'
            % weight by inverse standard deviation
            nTV = inputWeightsIn{2};
            nTV = nTV(:)';
            assert(length(nTV) == length(intensityV), 'mismatched length inputs');
            zIx = iswithintol(fractCorr, 0, eps*1e3);
            oneIx = iswithintol(fractCorr, 1, eps*1e3);
            tempC = fractCorr;
            tempC(zIx) = 0.5./nTV(zIx);  % mle is between 0 and 1/nTrials
            tempC(oneIx) = (nTV(oneIx)-0.5)./nTV(oneIx);
            inputWeights = sqrt(nTV ./ (tempC .* (1-tempC)));  
            
            if any(inputWeights > 1e10), keyboard; end
            debugPlot = false;
            if debugPlot
                oldFH = gcf;
                figure(100); clf;
                hold on;
                plot(inputWeights./sum(inputWeights), 'r');
                plot(nTV./sum(nTV), 'b');
            end
            outS.nTrials = nTV;
        otherwise
            error('Invalid inputWeights param');
    end
end    

%% some setup computations
midVal = mean(range(intensityV));

startingVals = [midVal, 1, fractCorr(end), fractCorr(1)];

optSet = optimset('Display', 'none');
inputWeights = inputWeights ./ sum(inputWeights);
outS.fitWeights = inputWeights;
if clampAtZero == 1
    modelFun = @(p,xs) (p(3) .* (1 - exp( -(xs./p(1)) .^ p(2))));    
    optFun = @(p) (modelFun(p,intensityV) - fractCorr).*inputWeights;
    startingVals = [midVal, 1, fractCorr(end)];

    [coefEsts resnorm residual exitflag] ...
        = lsqnonlin(optFun, startingVals, ...
                    [0 0 0], ...
                    [Inf Inf 1], ...
                    optSet);

    coefEsts(4) = 0;
elseif use50Thresh == 1
     modelFun = @(p,xs) (0.5+ (p(3)-0.5) .* (1 - exp( -(xs./p(1)) .^ p(2))));    
    optFun = @(p) (modelFun(p,intensityV) - fractCorr).*inputWeights;
    startingVals = [midVal, 1, fractCorr(end)];

    [coefEsts resnorm residual exitflag] ...
        = lsqnonlin(optFun, startingVals, ...
                    [0 0 0], ...
                    [Inf Inf 1], ...
                    optSet);

    coefEsts(4) = 0.5;
else
    modelFun = @(p,xs) (p(4) + (p(3)-p(4)) .* (1 - exp( -(xs./p(1)) .^ p(2))));
    optFun = @(p) (modelFun(p,intensityV) - fractCorr) .* inputWeights;
    startingVals = [midVal, 1, fractCorr(end), fractCorr(1)];
    
    [coefEsts resnorm residual exitflag] ...
        = lsqnonlin(optFun, startingVals, ...
                    [0 0 0 0], ...
                    [Inf Inf 1 1], ...
                    optSet);
end
    
if exitflag < 1
    disp('error in optimization');
    disp(exitflag)
end
                 
%% extract 50% threshold
p = coefEsts;
if use50Thresh
    oThresh = p(1) .* log(2).^(1/p(2));
    oThreshFractCorr = (p(3)-p(4)) ./ 2 + p(4);
else
    % default
    oThresh = p(1);
    oThreshFractCorr = 0.6321 * (p(3)-p(4)) + p(4);
end

% output
outS.thresh = oThresh;
outS.threshY = oThreshFractCorr;
outS.coefEsts = coefEsts;
outS.modelFun = modelFun;
outS.slope = coefEsts(2);    
outS.intensityX = intensityV;
outS.fractCorrY = fractCorr;



                 

