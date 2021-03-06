import math
import numpy
import random
def calc_tStimProbAvgLeft():
	probList = getvar('ProbList')
	doLowFirst = getvar('doLowProbFirst')
	numProbs = len(probList)

	compTrials = getvar('tNTrialsCompleted')
	nTrialPerBlock = getvar('RandProbTrialperBlock')
	lastIndex = getvar('tRandProbIndexLast')
	preblock = getvar('RandProbPreBlockTrials')
	blockNum = getvar('tRandBlockNumber')

	if compTrials == 0:
		if 0.5 in probList:
			ind = probList.index(0.5)
			setvar('tStimProbAvgLeft', 0.5)
			setvar('tRandProbIndexLast', ind)
			setvar('inPreBlock', 0)
			setvar('tRandBlockNumber',1)
		else:
			if doLowFirst == 1:
				tempval = 0
			else:
				tempval = numProbs-1
			setvar('tRandProbIndexLast', tempval)
			setvar('tRandBlockNumber',1)
			if preblock > 0:
				setvar('inPreBlock', 1)
				if probList[tempval] < 0.5:
					setvar('tStimProbAvgLeft', 0)
				if probList[tempval] > 0.5:
					setvar('tStimProbAvgLeft', 1)
			elif preblock == 0:
				setvar('inPreBlock', 0)
				setvar('tStimProbAvgLeft', probList[tempval])
	elif compTrials%nTrialPerBlock == 0:
		tempval = lastIndex
		blockNum = blockNum + 1
		setvar('tRandBlockNumber', blockNum)
		if doLowFirst == 0:
			if tempval == 1:
				if blockNum % (numProbs+1) == 0:
					tempval = 0;
					setvar('tRandProbIndexLast', tempval)
				else:
					tempval = numProbs-1;
					setvar('tRandProbIndexLast', tempval)
			else:
				tempval = 1
				setvar('tRandProbIndexLast', tempval)
		else:
			if tempval == 1:
				if blockNum % (numProbs+1) == 0:
					tempval = numProbs-1;
					setvar('tRandProbIndexLast', tempval)
				else:
					tempval = 0;
					setvar('tRandProbIndexLast', tempval)
			else:
				tempval = 1
				setvar('tRandProbIndexLast', tempval)	
		if preblock > 0:
			setvar('inPreBlock', 1)
			if probList[tempval] < 0.5:
				setvar('tStimProbAvgLeft', 0)
			if probList[tempval] > 0.5:
				setvar('tStimProbAvgLeft', 1)
			if probList[tempval] == 0.5:
				setvar('tStimProbAvgLeft', probList[tempval])
		elif preblock == 0:
			setvar('inPreBlock', 0)
			setvar('tStimProbAvgLeft', probList[tempval])
	elif preblock>0:
		if compTrials%nTrialPerBlock < preblock+1:
			setvar('inPreBlock', 1)
			setvar('tRandProbIndexLast', lastIndex)
			if probList[lastIndex] < 0.5:
				setvar('tStimProbAvgLeft', 0)
			if probList[lastIndex] > 0.5:
				setvar('tStimProbAvgLeft', 1)
			if probList[lastIndex] == 0.5:
				setvar('tStimProbAvgLeft', probList[lastIndex])
		if compTrials%nTrialPerBlock > preblock:
			setvar('inPreBlock', 0)	
			setvar('tStimProbAvgLeft', probList[lastIndex])
			setvar('tRandProbIndexLast', lastIndex)
	elif preblock == 0:
		setvar('inPreBlock', 0)		
		if compTrials%nTrialPerBlock > 0:
			setvar('tStimProbAvgLeft', probList[lastIndex])
			setvar('tRandProbIndexLast', lastIndex)


