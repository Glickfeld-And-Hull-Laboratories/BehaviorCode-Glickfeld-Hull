from matplotlib.pyplot import *
from numpy import *

import matplotlib.pyplot as plt
import numpy as np
import scipy.misc
import scipy.io
import scipy.signal.signaltools as s
import pdb
import os, os.path
import sys


outDName = '/tmp/measuredDpc.mat';
outD = scipy.io.loadmat(outDName, squeeze_me=True);

(nDpc, nLum) = shape(outD['meanLuminancePct'])
(crap, dpcNs) = unique(outD['gratingDpc'], return_index=True)

figure()
for iD in dpcNs:
    xv =  25/(2.0**(linspace(0,7,8)/(1.24)))
    semilogx(xv, outD['meanLuminancePct'][iD,:], '.-')

    xlabel('contrast (%)')
    ylabel('mean luminance change (%)')

legH = legend(outD['gratingDpc'][dpcNs], loc='lower left')

title('luminance change for stimuli as of 110112')

