import math
import numpy

def add_to_reactTimesCorrect():
	
	tReactTimeMsCorrect = getvar('actualHoldTimeMs') - getvar('tTotalReqHoldTimeMs')
	reactTimesMsCorrect = getvar('reactTimesMsCorrect')
	reactTimesMsCorrect.append(tReactTimeMsCorrect)
	setvar('reactTimesMsCorrect', reactTimesMsCorrect)
	targetCorrect = getvar('targetCorrect')

	mouseSpeed = getvar('mouseSpeed')

	if mouseSpeed == 0 and tReactTimeMsCorrect > 150 and tReactTimeMsCorrect < 550:
		targetCorrect = targetCorrect + 1
	elif mouseSpeed == 1 and tReactTimeMsCorrect > 250 and tReactTimeMsCorrect < 650:
		targetCorrect = targetCorrect + 1

	setvar('targetCorrect', targetCorrect)
	reactToWindow = targetCorrect/float(getvar('tTrialsDoneSinceStart'))
	setvar('reactToWindow', reactToWindow)

def add_to_reactTimesTotal():

    tReactTimeMsTotal = getvar('actualHoldTimeMs') - getvar('tTotalReqHoldTimeMs')
    reactTimesMsTotal = getvar('reactTimesMsTotal')
    reactTimesMsTotal.append(tReactTimeMsTotal)
    
    setvar('reactTimesMsTotal', reactTimesMsTotal)

def find_median():

	reactTimesMs = getvar('reactTimesMsTotal')

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

	randReqHoldMaxMs = getvar('randReqHoldMaxMs')

	reactToWindow = getvar('reactToWindow')
	
	reactToWindow_last40 = getvar('reactToWindow_last40')

	tTrialsDoneSinceStart = getvar('tTrialsDoneSinceStart')

	if reactToWindow >= 0.3 and (getvar('achievedTierRandom') == 0):
		setvar('randReqHoldMaxMs',randReqHoldMaxMs + 200)

	elif reactToWindow >= 0.5 and (getvar('achievedTierRandom') == 1) and (getvar('achievedMaxRandom') == 0):
		setvar('randReqHoldMaxMs',randReqHoldMaxMs + 400)

	elif reactToWindow <= 0.3 and randReqHoldMaxMs >= 200:
		setvar('randReqHoldMaxMs',randReqHoldMaxMs)
	#	setvar('randReqHoldMaxMs',randReqHoldMaxMs - 200)

	elif (reactToWindow + reactToWindow_last40)/2 <= 0.4 and (tTrialsDoneSinceStart%80 == 0) and randReqHoldMaxMs >= 1200:
		setvar('randReqHoldMaxMs',randReqHoldMaxMs)
		#setvar('randReqHoldMaxMs',randReqHoldMaxMs - 400)

	if getvar('randReqHoldMaxMs') >= 1200:
		setvar('achievedTierRandom', 1)

	elif getvar('randReqHoldMaxMs') == 4000:
		setvar('achievedMaxRandom', 1)

	if tTrialsDoneSinceStart % 40 == 0:
		setvar('reactToWindow_last40', reactToWindow)


def updateReactTime():	

	if getvar('mouseHolder') == 1 and getvar('reactTimeMs') >= 2000:
		reactTimeMs = getvar('reactTimeMs')
		setvar('reactTimeMs', reactTimeMs - 1000)

	if (getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 0 and reactToWindow >= 0.6 and getvar('reactTimeMs') != 600:
		setvar('reactTimeMs', 600)

	elif(getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 1 and reactToWindow >= 0.6 and getvar('reactTimeMs') != 720:
		setvar('reactTimeMs', 720)

	elif (getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 0 and reactToWindow >= 0.45 and getvar('reactTimeMs') != 1000:
		setvar('reactTimeMs', 1000)

	elif (getvar('achievedMaxRandom') == 1) and getvar('mouseSpeed') == 1 and reactToWindow >= 0.45 and getvar('reactTimeMs') != 1200:
		setvar('reactTimeMs', 1200)
