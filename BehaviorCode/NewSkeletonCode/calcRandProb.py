import math
import numpy
import random
def calc_tStimProbAvgLeft():
	probList = getvar('ProbList')
	numProbs = len(probList)

	compTrials = getvar('tNTrialsCompleted')
	nTrialPerBlock = getvar('RandProbTrialperBlock')
	lastIndex = getvar('tRandProbIndexLast')

	if compTrials == 0:
		if 0.5 in probList:
			ind = probList.index(0.5)
			setvar('tStimProbAvgLeft', 0.5)
			setvar('tRandProbIndexLast', ind)
		else:
			tempval = random.randint(1,numProbs)-1
			setvar('tStimProbAvgLeft', probList[tempval])
			setvar('tRandProbIndexLast', tempval)
	elif compTrials%nTrialPerBlock > 0:
		print lastIndex
		setvar('tStimProbAvgLeft', probList[lastIndex])
	elif compTrials%nTrialPerBlock == 0:
		tempval = lastIndex
		while(tempval == lastIndex):
			tempval = random.randint(1,numProbs)-1
		setvar('tStimProbAvgLeft', probList[tempval])
		setvar('tRandProbIndexLast', tempval)


