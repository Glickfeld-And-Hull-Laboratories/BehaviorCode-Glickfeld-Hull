import math
import numpy

def calc_tStimProbAvgLeft():
	probList = getvar(ProbList)
	numProbs = len(probList)

	compTrials = getvar(tNTrialsCompleted)
	nTrialPerBlock = getvar(RandProbTrialperBlock)
	lastIndex = getvar(tRandProbIndexLast)

	if compTrials == 0
		if 0.5 in probList
			ind = probList.index('0.5')
			setval('tStimProbAvgLeft', 0.5)
			setval('tRandProbIndexLast', ind)
		else
			tempval = random.randrange(1,numProbs)
			setval('tStimProbAvgLeft', probList(tempval))
			setval('tRandProbIndexLast', probList(tempval))
	elif compTrials%nTrialPerBlock > 0
		setval('tStimProbAvgLeft', probList(lastIndex))
	elif compTrials%nTrialPerBlock == 0
		tempval = lastIndex
		while(tempval == lastIndex)
		tempval = random.randrange(1,numProbs)
		setval('tStimProbAvgLeft', probList(tempval))
		setval('tRandProbIndexLast', probList(tempval))


