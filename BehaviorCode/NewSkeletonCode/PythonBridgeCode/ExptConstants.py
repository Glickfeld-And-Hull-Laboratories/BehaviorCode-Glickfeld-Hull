strobedWords = { 
    'START' : 170, 
    'END' : 85, 
    'TimestampStart' : 200, 
    'TimestampEnd' : 201, 
    'LeverPressStart' : 3, 
    'VisStimOn' : 4, 
    'LaserStimOn' : 5, 
    'ItiStart' : 6, 
    'StartTrialWaitForPress' : 7, 
    'LeverReleaseEarly' : 8, 
    'LeverReleaseCorrect' : 9, 
    'Reward' : 10, 
    'LeverNoReleaseMiss' : 11, 
    'LeverSolenoidOn' : 12, 
    'LeverSolenoidOff' : 13, 
    'VisStimOff' : 14 }

experimentXmlTrialIds = { 
    'HoldAndDetectConstant8' :  8, 
    'HoldAndDetectConstant5' :  10, 
    'DirTuningMapping' :  11, 
    'ChRMappingTwelve' :  12, 
    'ChRMappingTen' :  13, 
    'DualPress' :  14,
    'Lego' : 15, 
    'VisStimRet' : 90}

trialOutcomes = [ 
    'success', 
    'early', 
    'ignore', 
    'failure', 
    'incorrect', 
    'dualrelease' ];

reverseExperimentXmlTrialIds = dict(zip(experimentXmlTrialIds.values(), experimentXmlTrialIds.keys()))

