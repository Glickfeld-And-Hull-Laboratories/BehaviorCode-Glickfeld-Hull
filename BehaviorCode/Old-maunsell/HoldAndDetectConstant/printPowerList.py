import numpy as np

tMax = 10
nSteps = 8
stepsPerOct = 2

powerListMw = tMax / (2.0 ** (np.arange(0,nSteps, dtype='float')/stepsPerOct))
print powerListMw
