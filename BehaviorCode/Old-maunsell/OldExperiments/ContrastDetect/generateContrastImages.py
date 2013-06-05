#!/Library/Frameworks/Python.framework/Versions/Current/bin/python2.6
#  above line ensures we use the EPD python distribution
#
#  Note contrast is specified as peak-to-peak; contrast of 100% is a sinusoid
#  which goes from 0...255 in pixel intensity; contrast of 50% is 63-191, etc.
#
#  Edit the parameter file (default):
#     parameterFile = '~/Library/Application Support/MWorks/ContrastDetectGeneratedImages/imageParameters.py'
#  to set up the image parameters, then run this script
#
#  created: histed 2010/7/16


####### Don't change these
import os

parameterFile = 'imageParameters.py'  # default in outDir
outDir = os.path.join(os.getenv('HOME'), 'Library/Application Support/MWorks/ContrastDetectGeneratedImages')
#outDir = './generated-images'
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

# get parameters
sys.path.append(outDir)
import imageParameters as p

# generate contrast vector
powRange = np.arange(0, p.nContrasts*1.0, 1) / p.stepsPerOctave;
contrastPct = p.maxContrastPct / np.power(2, powRange);

# status
print 'Generating images:'
print '  %d images ' % p.nContrasts
print '  contrast levels (percents): %s ' % np.array_str(contrastPct, precision=2)


################
# subfuncs
def safe_mkdir(tPath):
    try:
        os.makedirs(tPath)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

################

# generate images
nImages = p.nContrasts
monitorFullWidthDeg = np.arctan2(p.monitorWidthCm, p.monitorDistanceCm) / np.pi * 180
print "Size of monitor %5.1f deg - Note this is wrong!! 110112 MH" % (monitorFullWidthDeg)
#monitorFullHeightDeg = np.arctan2(monitorHeightCm, monitorDistanceCm) / np.pi * 180
pixPerDegX = (p.monitorWidthPix/monitorFullWidthDeg) #, monitorHeightPix/monitorFullHeightDeg)

for iI in range(0, nImages):
    tContrastPct = contrastPct[iI]

    # make an image in floats, 0-1
    xyValsPix = np.mgrid[0:2*p.monitorWidthPix,0:2*p.monitorHeightPix]

    # below we use only monitor X
    xValsRad = xyValsPix[0] / pixPerDegX / 180 * np.pi * 2 # *2 because we double size above

    # create a vertically-oriented grating
    pixVals = np.sin( 1.0/(float(p.gratingDPC)/180*np.pi) * 2*np.pi * xValsRad)

    # now pixVals ranges from -1 to 1
    # rescale to 0..256, using correct contrast
    intPixVals = np.asarray( (pixVals*128) * tContrastPct/100 + 127, dtype='uint8') 

    # rotate to proper direction
    rotPixVals = scipy.misc.imrotate(intPixVals, p.gratingOriDeg-90)  # rotate back to 0 deg first
    (nX, nY) = np.shape(rotPixVals)
    rotPixVals = rotPixVals[p.monitorWidthPix/2:p.monitorWidthPix/2+p.monitorWidthPix,
                            p.monitorHeightPix/2:p.monitorHeightPix/2+p.monitorHeightPix]
    
    # rescale
    rotPixVals = np.transpose(rotPixVals)
    rgbPixVals = np.dstack( (rotPixVals,rotPixVals,rotPixVals) )

    # save image
    im = scipy.misc.toimage(rgbPixVals)
    outDir = os.path.expanduser(outDir)
    outName = os.path.join(outDir, 'generated-%d.png' % iI)
    safe_mkdir(outDir)
    scipy.misc.imsave(outName, im)

    print 'Done, image %d of %d, contrast %4.1f%%' % (iI+1, nImages, tContrastPct)

# save image data to a mat file
datOutName = os.path.join(outDir, dataOutName)
outD = {}
outD['runTime'] = datetime.isoformat(datetime.now(), ' ')
outD['sourceProgramName'] = sys.argv[0]
outD['nImages'] = nImages
outD['contrastLevelsPct'] = contrastPct
for tn in ['maxContrastPct','nContrasts','stepsPerOctave','gratingType','gratingOriDeg',
           'gratingDPC', 'monitorDistanceCm', 'pixPerDegX', 'monitorFullWidthDeg', 
           'monitorWidthCm']:
    try: 
        outD[tn] = eval(tn)
    except:
        outD[tn] = eval('p.'+tn)
    
scipy.io.savemat(datOutName, outD, oned_as='column')
print 'Saved matlab info data to "%s"' % datOutName



