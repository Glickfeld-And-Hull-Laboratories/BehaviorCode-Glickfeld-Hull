import math
import numpy
import random
def calc_tGratingDirectionDeg():
	dirList = getvar('gratingDirections')
	numDirs = len(dirList)

	dDirList = getvar('dGratingDirections')
	numdDirs = len(dDirList)

	tempval1 = random.randint(1,numDirs) - 1
	setvar('tGratingDirectionDeg', dirList[tempval1])

	if numdDirs > 1:
		tempval2 = random.randint(1,numdDirs) - 1
		setvar('dGratingDirectionDeg', dDirList[tempval2])
	else:
		setvar('dGratingDirectionDeg', dDirList[0])
	

	rewList = getvar('rewardStim')

	if dirList[tempval1] in rewList:
		setvar('tRewardTrial', 1)
	else:
		setvar('tRewardTrial', 0)


