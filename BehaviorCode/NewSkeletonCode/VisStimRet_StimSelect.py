import math
import numpy
import random

usedList = getvar('usedList')

nStimCond = getvar('nStimCond')
stimList = list(range(0,nStimCond*5))

unUsedList = list(set(usedList).symmetric_difference(set(stimList)))
if len(unUsedList)==0:
	unUsedList = stimList
	usedList = []

doRand = getvar('doRand')
if doRand == 1:
	selectStim = random.choice(unUsedList)
else:
	selectStim = unUsedList(0)

usedList.append(selectStim)

setvar('usedList',usedList)
setvar('tStimulusNumber',selectStim)

