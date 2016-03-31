import math
import numpy

def add_to_reactTimesCorrect():

    tReactTimeMsCorrect = getvar('actualHoldTimeMs') - getvar('tTotalReqHoldTimeMs')
    
    reactTimesMsCorrect.append(tReactTimeMsCorrect)

    setvar('reactTimesMsCorrect', reactTimesMsCorrect)

    if tReactTimeMsCorrect > 150 and tReactTimeMsCorrect < 550:
    	fastReact = getvar('fastReact') + 1
    	setvar('fastReact', fastReact)

    elif tReactTimeMsCorrect > 250 and tReactTimeMsCorrect < 650:
    	slowReact = getvar('slowReact') + 1
    	setvar('slowReact', slowReact) 

def add_to_reactTimesTotal():

    tReactTimeMsTotal = getvar('actualHoldTimeMs') - getvar('tTotalReqHoldTimeMs')
    
    reactTimesMsTotal.append(tReactTimeMsTotal)
    
    setvar('reactTimesMsTotal', reactTimesMsTotal)

def find_median():

	reactTimesMs = getvar('reactTimesMs')

	medianReactTimesMs = numpy.median(reactTimesMs)

	setvar('medianReactTimesMs', medianReactTimesMs)

	if medianReactTimesMs < 400:
		setvar('mouseSpeed', 0)
	else:
		setvar('mouseSpeed', 1)

	if medianReactTimesMs < 750:
		setvar('mouseHolder', 0)
	else:
		setvar('mouseHolder', 1)

def get_targetCorrect():

	targetCorrect = getvar('fastReact')/float(getvar('slowReact'))
	randReqHoldMaxMs = getvar('randReqHoldMaxMs')
	targetCorrect_last40 = getvar('targetCorrect_last40')
	
	if targetCorrect >= 0.3 and (getvar('achievedTierRandom') == 0):
		setvar('randReqHoldMaxMs',randReqHoldMaxMs + 200)

	elif targetCorrect <= 0.3:
		setvar('randReqHoldMaxMs',randReqHoldMaxMs - 200)

	elif targetCorrect >= 0.55 and (getvar('achievedMaxRandom') == 0):
		setvar('randReqHoldMaxMs',randReqHoldMaxMs + 400)

	elif (targetCorrect + targetCorrect_last40)/2 <= 0.4 and (getvar('tTrialsDoneSinceStart')%80 == 0):
		setvar('randReqHoldMaxMs',randReqHoldMaxMs - 200)

	if getvar('randReqHoldMaxMs') == 1200:
		setvar('achievedTierRandom', 1)

	elif getvar('randReqHoldMaxMs') == 4000:
		setvar('achievedMaxRandom', 1)

	targetCorrect_last40 = targetCorrect
	setvar('targetCorrect_last40', targetCorrect_last40)
	setvar('fastReact', 0)
	setvar('slowReact', 0)

def updateReactTime():

	if getvar('mouseHolder') == 1:
		reacTimeMs = getvar('reacTimeMs')
		setvar('reacTimeMs', reacTimeMs - 1000)

	if (getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 0 and getvar('reacTimeMs') != 600:
		setvar('reacTimeMs', 600)

	elif (getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 1 and getvar('reacTimeMs') != 720:
		setvar('reacTimeMs', 720)
