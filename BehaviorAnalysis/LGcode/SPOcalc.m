function [SPO, out] = SPOcalc(maxStep,minStep,nStep)
% SPOcalc gives the number of steps per octave (SPO) and all intermediate 
% steps (out) when given a min, max and number of steps desired
% [SPO, out] = SPOcalc(maxStep,minStep,nStep)

SPO = (nStep-1)./(log2(maxStep./minStep));
out = zeros(1,nStep);
for iStep = 1:nStep
    out(1,iStep) = maxStep./(2.^((iStep-1)/SPO));
end