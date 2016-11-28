import math
import numpy

def calc_tStimProbAvgLeft():
	probList = getvar('ProbList')
	numProbs = len(probList)

	compTrials = getvar('tNTrialsCompleted')
	nTrialPerBlock = getvar('RandProbTrialperBlock')
	lastIndex = getvar('tRandProbIndexLast')

	if compTrials == 0
		if 0.5 in probList
			ind = probList.index('0.5')
			setvar('tStimProbAvgLeft', 0.5)
			setvar('tRandProbIndexLast', ind)
		else
			tempval = random.randrange(1,numProbs)
			setvar('tStimProbAvgLeft', probList(tempval))
			setvar('tRandProbIndexLast', probList(tempval))
	elif compTrials%nTrialPerBlock > 0
		setvar('tStimProbAvgLeft', probList(lastIndex))
	elif compTrials%nTrialPerBlock == 0
		tempval = lastIndex
		while(tempval == lastIndex)
			tempval = random.randrange(1,numProbs)
		setvar('tStimProbAvgLeft', probList(tempval))
		setvar('tRandProbIndexLast', probList(tempval))


