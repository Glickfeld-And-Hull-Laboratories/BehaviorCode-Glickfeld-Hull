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


if len(sys.argv) < 3: doExport = False

outDir = os.path.join(os.getenv('HOME'), 'Library/Application Support/MWorks/ContrastDetectGeneratedImages')

def doCompute(doPlot=False, fileNum=7): 
    b1 = scipy.misc.imread(os.path.join(outDir,'generated-%d.png' % fileNum))

    r1Orig = np.array(b1[0,:,0], dtype='f8')
    
    
    # reinterp to a lot
    nInterp = 10
    nPix = len(r1Orig)
    r1Interp = scipy.interp(linspace(0,nPix,nPix*nInterp), range(nPix), r1Orig)
    print len(r1Interp)

    
    d1 = diff(r1Interp)
    d2 = diff(d1)


    # find zero crossings
    zcIx = diff(sign(r1Interp-127))>=1
    zcNs = where(zcIx)[0]
    tIx = logical_and(zcIx,d1> 1.0/nInterp/2)
    tNs = nonzero(tIx)[0]

    # reduce by one so it matches the diff length
    #r1 = r1Interp[:-1]  


    if doPlot:
        figure()
        plot(r1Interp[:-1])
        xv = arange(nPix*nInterp)[:-1]
        plot(xv[tIx], r1Interp[:-1][tIx], 'rx')
        plot((nPix*nInterp)/2 * array((1,1)), ylim(), 'g')


    screenPix = 1280
    screenDeg = 72

    periodPix = mean(diff(tNs))/nInterp
    if len(tNs) < 2 or ptp(diff(tNs)) > periodPix*nInterp/5:
        print '**error calculating zero crossings: numbers may be wrong!'
        isGood = False

        #import pdb; pdb.set_trace()

        figure()
        plot(r1Interp[:-1])
        xv = arange(nPix*nInterp)[:-1]
        plot(xv[tIx], r1Interp[:-1][tIx], 'rx')
        plot((nPix*nInterp)/2 * array((1,1)), ylim(), 'g')

    else:
        isGood = True

    dpc = screenDeg / (screenPix/periodPix)

    # compute mean luminance change
    mL = mean(r1Orig.flat)  # not interp
    pctChange = (mL-127.0)/127.0*100

    # summary
    print('Image %2d - Measured DPC: %5.1f, mean luminance change from 127: %5.2f%%'
          % (fileNum, dpc, pctChange))

    return (dpc, pctChange, isGood, nPix)

##################

#doCompute(True, 4)

if len(sys.argv) > 1:
    doPlot = sys.argv[1]
else:
    doPlot = False
    
allDpc = zeros(8);
allPct = zeros(8);
allGood = zeros(8);
for iA in range(8):
    (allDpc[iA], allPct[iA], allGood[iA], monitorWidthPix) = doCompute(doPlot, iA)

# some computations
estDpc = mean(allDpc[:5])

# get data and save it
paramD = scipy.io.loadmat(os.path.join(outDir, 'info.mat'))


# 
outDName = '/tmp/measuredDpc.mat';
if not os.path.isfile(outDName):
    # init it
    outD = {}
    outD['gratingDpc'] = empty((0),dtype='f8')
    outD['monitorWidthPix'] = empty((0),dtype='f8')
    outD['monitorWidthCm'] = empty((0),dtype='f8')
    outD['measuredDpc'] = array([])
    outD['meanLuminancePct'] = zeros([0,8])

else:
    # load it
    outD = scipy.io.loadmat(outDName, squeeze_me=True, struct_as_record=False);

# add data
outD['gratingDpc'] = hstack((outD['gratingDpc'], paramD['gratingDPC'][0]))
outD['monitorWidthPix'] = hstack((outD['monitorWidthPix'], asarray(monitorWidthPix, dtype='f')))
outD['monitorWidthCm'] = hstack((outD['monitorWidthCm'], paramD['monitorWidthCm'][0]))
outD['measuredDpc'] = hstack((outD['measuredDpc'], asarray(estDpc)))
outD['meanLuminancePct'] = vstack((outD['meanLuminancePct'], asarray(allPct.reshape((1,8)))))

# save it back
if doExport:
    os.system('mv %s %s-bak 2>&1 > /dev/null' % (outDName, outDName))
    scipy.io.savemat(outDName, outD, oned_as='column')
