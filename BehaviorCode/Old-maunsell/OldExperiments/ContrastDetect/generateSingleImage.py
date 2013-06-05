#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
#  above line ensures we use the EPD python distribution
#
#  created: histed 2010/01/16

# todo: mean should be 127, should be odd-symmetric


####### Don't change these
import os

#parameterFile = 'imageParameters.py'  # default in outDir
outDir = '/tmp/'
dataOutName = 'info.mat'  # in above directory

################


# imports
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
import scipy.misc
import PIL
import errno, sys
from datetime import datetime
import scipy.io
from numpy import *
from matplotlib.pyplot import *

# get parameters
pixPerCycle = 4   # width 1280 gives 320 c/screen; 72 deg screen gives 4.5 cpd
#screenWidthPix = 1280
#screenHeightPix = 1024
screenWidthPix = 2048#2048
screenHeightPix = 2048

################
# subfuncs

################

# generate images
#xs = linspace(-pi, pi, screenWidthPix)
nCyclesHalf = screenWidthPix/pixPerCycle/2
xs1 = r_[-(screenWidthPix/2):screenWidthPix/2]
xs = xs1/float(pixPerCycle)*2*pi

xyN = np.mgrid[0:screenHeightPix,0:screenWidthPix]

pixVals = np.cos(xs[xyN[1,:,:]])

pixF = array(np.round(pixVals*127+127), dtype='uint8')
#pixF = array(np.round(pixVals*126+127), dtype='uint8')

rgbPixVals = np.dstack((pixF,pixF,pixF))

# display test
clf()
plt.imshow(pixF, cmap=cm.gray)

# save image
im = scipy.misc.toimage(rgbPixVals)
outName = os.path.join(outDir, 'generated-single.png')
scipy.misc.imsave(outName, im)

mI = mean(np.ravel(im))
print 'Done, image has mean %4.2f' % mI

## # save image data to a mat file
## datOutName = os.path.join(outDir, dataOutName)
## outD = {}
## outD['runTime'] = datetime.isoformat(datetime.now(), ' ')
## outD['sourceProgramName'] = sys.argv[0]
## outD['nImages'] = nImages
## outD['contrastLevelsPct'] = contrastPct
## for tn in ['maxContrastPct','nContrasts','stepsPerOctave','gratingType','gratingOriDeg',
##            'gratingDPC', 'screenDistanceCm', 'pixPerDegX', 'screenFullWidthDeg', 
##            'screenWidthCm']:
##     try: 
##         outD[tn] = eval(tn)
##     except:
##         outD[tn] = eval('p.'+tn)
    
## scipy.io.savemat(datOutName, outD, oned_as='column')
## print 'Saved matlab info data to "%s"' % datOutName



