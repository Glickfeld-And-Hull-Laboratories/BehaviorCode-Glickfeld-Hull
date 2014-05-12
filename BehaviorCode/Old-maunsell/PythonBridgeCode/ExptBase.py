# ExptBaseMH.py:
# Shared code to handle experiment duties
#       Talks to Blackrock over serial; to laser servers over network
#       Chooses constants/functionality based on machine hostname
#       Registers callbacks for most variables and maintains their current state
#
# MH 130102: broke out from HoldAndDetectLaserPython.py

import sys
import socket
import sys
import cPickle as p
import struct
import time
import numpy as np
import math
import threading
import socket
import string
import copy
import gzip
from cStringIO import StringIO
from copy import deepcopy
import ExptConstants as cs

sys.path.append('/Library/Application Support/MWorks/Scripting/Python')

from mworks.conduit import IPCClientConduit

def chop(x, sig=2):
    nSig = sig-(np.floor(np.log10(x)))-1
    if x == 0:
        return 0
    if np.size(x) == 1:
        x = [x]
        nSig = [nSig]
    chopped = np.array([np.around(x0,int(n0)) for (x0,n0) in zip(x,nSig)])
    return chopped


################################################################
## set machine-specific constants

thisHostname = socket.gethostname()
thisHostname = string.replace(thisHostname, '.dhe.duke.edu', '')
if thisHostname == "hullglicks-mini":
    PORT = 9991
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouse2":
    PORT = 9991
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouse3":
    PORT = 9992
    PORT2 = 9993  # use same as rig 4
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouse4":
    PORT = 9993
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouseTest":
    PORT = None
    doTransmitCodesToBlackrock = False
else:
    raise RuntimeError, 'Found unknown hostname: %s' % thisHostname
    
HOST = "152.16.224.214"
laserParamCode = 'sendLaserParams'
lastSendTime = 0
nCodesSent = 0

if doTransmitCodesToBlackrock:
    import json
    import serial

################################################################
## state-retaining class

class stateC:
    lock = None
    nCodes = None
    nStrobedTotal = 0
    eventDict = None
    experimentXmlTrialId = None
    blankEventDict = { 'names':[], 'values':[], 'ts':[] }
    variableCurrValues = {}
    variableLastTs = {}
    variableThisTrialValues = {}
    variableThisTrialTs = {}
    clientAccessLock = None
    codecCopy = None
    requiredExptId = None # set on initialization, must never change

    def __init__(self, exptIdStr):
        self.lock = threading.Lock()
        self.clientAccessLock = threading.Lock()

        self.eventDict = deepcopy(self.blankEventDict)
        self.requiredExptId = cs.experimentXmlTrialIds[exptIdStr]  # checked in cbParseVariables
        
        
    def getEventDict(self):
        """return it and clear it"""
        self.lock.acquire()
        edCopy = self.eventDict.copy()
        self.eventDict = deepcopy(self.blankEventDict)
        self.lock.release()
        return edCopy

    def appendToEventDict(self, name, data, time):
        self.lock.acquire()
        self.eventDict['names'].append(name)
        self.eventDict['values'].append(data)
        self.eventDict['ts'].append(time)
        self.lock.release()

    def saveVariableValue(self, name, evts):
        self.lock.acquire() # needed only during writes to data structures
        self.variableCurrValues[name] = evts.data
        self.variableLastTs[name] = evts.time
        self.variableThisTrialValues[name] = evts.data
        self.variableThisTrialTs[name] = evts.time
        self.lock.release()

    def getVariableValues(self):
        return copy.deepcopy(self.variableCurrValues)

    def resetThisTrialData(self):
        self.lock.acquire()
        self.variableThisTrialValues = {}
        self.variableThisTrialTs = {}
        self.lock.release()


################################################################
## functions


def computeLaserPowerHADC8():
    """pass in variableCurrValues,
    returns a dictionary that can be passed over the wire to the laser stim machine"""
    inD = state.variableCurrValues
    laserTrainRandomType = None
    
    if (inD['laserTrainRandomType'] == 0 or inD['laserTrainRandomType'] == 0.0
        or len(inD['laserTrainRandomType']) == 0):
        laserTrainRandomType = None
    else:
        laserTrainRandomType = inD['laserTrainRandomType']   # string

    if laserTrainRandomType == 'uniformRandom':
        nPulses = inD['laserTrainRandomNPulses']
    else:
        nPulses = None

        
    outD = { 'pulsePowerMw': float(inD['tLaserPowerMw']),
             'offPowerMw': float(inD['laserOffPowerMw']),
             'baselinePowerMw': float(inD['tLaserBaselinePowerMw']),
             'baselinePostStimTimeMs': float(inD['laserBaselinePostStimTimeMs']),
             'baselineDelayTimeS': 0.0, #don't turn off in ITI: float(inD['tItiWaitTimeMs']-500)/1000, 
             'doLinearRamp': int(inD['tLaserDoLinearRamp']),
             'rampLengthMs': float(inD['tLaserRampLengthMs']),
             'doPulseTrain': int(inD['tLaserDoPulseTrain']),
             'pulseLengthMs': float(inD['laserPulseLengthMs']),
             'pulsePeriodMs': float(inD['laserPulsePeriodMs']),  
             'trainRandomType': laserTrainRandomType, # string or None
             'trainNPulses': nPulses, # only used if random
             'trainLengthMs': float(inD['laserTrainLengthMs']),
             'trainRampUpMs' : float(inD['laserTransitionRampUpDownMs']),
             'trainRampDownMs' : float(inD['laserTransitionRampUpDownMs']),
             'trainRampUpDownDoExp': float(inD['laserTransitionDoExpRamp']),
             'rampRampDownMs' : float(inD['laserTransitionRampUpDownMs']),             
             'rampRampUpDownDoExp' : float(inD['laserRampDoExpRamp']),
             'rampExtraConstantLengthMs' : float(inD['laserRampExtraConstantLengthMs']),
             }
    assert(outD['pulsePowerMw'] >= 0)

    return outD

def computeLaserPowerTrialLaser(doDirTuningMapping=False):
    """as above"""
    inD = state.variableCurrValues

    if not doDirTuningMapping:  # normal HandDC8 - use time of trial
        #totalTrLen = inD['fixedReqHoldTimeMs']+inD['randReqHoldMaxMs']+inD['reactTimeMs']
        thisTrLen = inD['tTotalReqHoldTimeMs']+inD['reactTimeMs']  # minimize the time laser is on
    else:  # dir tuning - use stimulus on time
        thisTrLen = (inD['stimulusOnTimeMs'] + inD['postStimulusTimeMs'] \
                     + (inD['itiTimeMs']-inD['laserOnTimeFromStartOfItiMs']))
        
    if inD['tTrialLaserOnTimeMs'] == 0 and inD['tTrialLaserOffTimeMs'] == 0:
        pulseLen = thisTrLen
        pulsePeriod = 0
    else:
        pulseLen = inD['tTrialLaserOnTimeMs']
        pulsePeriod = inD['tTrialLaserOnTimeMs'] + inD['tTrialLaserOffTimeMs']

    outD = {
        'pulsePowerMw': float(inD['tTrialLaserPowerMw']),
        'offPowerMw': float(inD['laserOffPowerMw']),
        'doLinearRamp': int(0),
        'doPulseTrain': int(1),
        'pulseLengthMs': float(pulseLen), 
        'pulsePeriodMs': float(pulsePeriod),
        'trainRandomType': None, # not supported unless random
        'trainNPulses': None,  # not supported unless random
        'trainLengthMs': float(thisTrLen), 
        'trainRampUpMs' : float(inD['laserTransitionRampUpDownMs']),
        'trainRampDownMs' : float(inD['laserTransitionRampUpDownMs']),
        'trainRampUpDownDoExp': float(inD['laserTransitionDoExpRamp']),
        }
    return outD
        
# change to do dir mapping here
def computeLaserPowerChRMap():
    """as above"""
    inD = state.variableCurrValues

    if inD['tRampTrain'] == 10:
        # ramp
        doLinearRamp = 1
        doPulseTrain = 0
        trainLengthMs = 0
    elif inD['tRampTrain'] == 11:
        #train
        doLinearRamp = 0
        doPulseTrain = 1
        tPulseLengthMs = inD['tTrainPulseLengthMs']+2*inD['laserTransitionRampUpDownMs']
        trainLengthMs = (inD['tTrainNPulses']-1) * inD['tTrainPeriodMs'] + tPulseLengthMs
        
    outD = { 'pulsePowerMw': float(inD['tPeakPowerMw']),
             'offPowerMw': float(inD['laserOffPowerMw']),
             'doLinearRamp': int(doLinearRamp), 
             'rampLengthMs': float(inD['tRampLengthMs']),
             'doPulseTrain': int(doPulseTrain), 
             'pulseLengthMs': float(inD['tTrainPulseLengthMs']),
             'pulsePeriodMs': float(inD['tTrainPeriodMs']),
             'trainRandomType': None, # not supported unless random
             'trainNPulses': None,  # not supported unless random
             'trainLengthMs': float(trainLengthMs), 
             'trainRampUpMs' : float(inD['laserTransitionRampUpDownMs']),
             'trainRampDownMs' : float(inD['laserTransitionRampUpDownMs']),
             'trainRampUpDownDoExp': float(inD['laserTransitionDoExpRamp']),
             'rampRampDownMs' : float(inD['laserTransitionRampUpDownMs']),
             'rampRampUpDownDoExp' : float(inD['laserTransitionDoExpRamp']),
             'rampExtraConstantLengthMs' : float(inD['tRampExtraConstantLengthMs']),    
             }
    assert(outD['pulsePowerMw'] >= 0)

    return outD
    

def cbSendDigitalWordSerial(evt):
    state.appendToEventDict('dW', evt.data, evt.time)

    ### debug
    ## if evt.data == 4:
    ##     print 'word', evt.data
    ##     print evt.time


    state.nStrobedTotal += 1
    if state.nStrobedTotal % 100 == 0:
        print("Received %d strobed codes total" % state.nStrobedTotal)
    sys.stdout.flush()

def cbSendLaserParams(evts):
    #print evts.data; sys.stdout.flush()
    if evts.data == 1:
        # do only on true
        if (state.experimentXmlTrialId == cs.experimentXmlTrialIds['HoldAndDetectConstant5']
            or state.experimentXmlTrialId == cs.experimentXmlTrialIds['HoldAndDetectConstant8']):

            if ((state.variableCurrValues['doBlock2'] == 1 and state.variableCurrValues['block2DoTrialLaser'] == 1)
                or (state.variableCurrValues['trialLaserPowerMw'] > 0)):
                laserOutDict = computeLaserPowerTrialLaser()
            else:
                laserOutDict = computeLaserPowerHADC8()
            sendLaserParams(laserOutDict)

        elif (state.experimentXmlTrialId == cs.experimentXmlTrialIds['DirTuningMapping']):
              if state.variableCurrValues['tTrialLaserPowerMw'] > 0:
                # else: if power is zero, means no laser this trial, just skip
                laserOutDict = computeLaserPowerTrialLaser(doDirTuningMapping=True)
                sendLaserParams(laserOutDict)

        elif (state.experimentXmlTrialId == cs.experimentXmlTrialIds['ChRMappingTwelve'])\
             or (state.experimentXmlTrialId == cs.experimentXmlTrialIds['ChRMappingTen']):
            #print variableCurrValues
            laserOutDict = computeLaserPowerChRMap()
            sendLaserParams(laserOutDict)

        elif (state.experimentXmlTrialId == cs.experimentXmlTrialIds['DualPress']):
            print "DualPress laser not implemented yet; not sending to server"
            # todo: MH 130102

        else:
            raise RuntimeError, 'trial Id %d not matched' % state.experimentXmlTrialId
        

def cbSendSerialParams(evts):
    """do actual serial transmission in this callback fcn"""
    global sp
    global client
    global lastSendTime
    global nCodesSent

    if evts.data == 1: # do only on true
        # check for error where we get two send events w/ same timestamp
        if lastSendTime == evts.time:
            raise RuntimeError, 'Two send events in a row found'
            #print evts.time
            #print evts.data
        else:
            lastSendTime = evts.time

        # take a snapshot of values to send
        serialOutDict = deepcopy(state.variableCurrValues)
        stateOutDict = state.getEventDict();

        # debug checks for variable presence
        if state.experimentXmlTrialId == cs.experimentXmlTrialIds['HoldAndDetectConstant5']:
            if not serialOutDict.has_key('tGratingContrast'):
                raise RuntimeError, 'Missing tGratingContrast'
        elif state.experimentXmlTrialId == cs.experimentXmlTrialIds['DirTuningMapping']:
            if not serialOutDict.has_key('tGratingDirectionDeg'):
                raise RuntimeError, 'Missing tGratingDirectionDeg'
        elif (state.experimentXmlTrialId == cs.experimentXmlTrialIds['ChRMappingTwelve']
              or state.experimentXmlTrialId == cs.experimentXmlTrialIds['ChRMappingTen']):
            chkList = ['tRampTrain', 'tTrainNPulses', 'tTrainPeriodMs', 'tTrainPulseLengthMs',
                       'juice', 'tDuringStimTimeMs', 'tGiveRewardThisTrial', 'tIsRewardGiven',
                       'tPostStimTimeMs', 'tRewardTimeFromItiOnMs', 'tPeakPowerMw' ]
            for c in chkList:
                if not serialOutDict.has_key(c):
                    raise RuntimeError, 'Missing variable: %s' % c
        
        # first summary data
        outStr = json.dumps(serialOutDict, separators=(',',':')) + '\n'
        # 2nd codec
        state.clientAccessLock.acquire()
        # was: client.codec
        outStr = outStr+json.dumps(state.codecCopy, separators=(',',':')) + '\n'
        state.clientAccessLock.release()
        # 3rd, event stream
        outStr = outStr + json.dumps(stateOutDict, separators=(',',':')) + '\n'

        # compress it
        # 120920: does not work, need a way to delimit start and end of each gzip stream
        doCompress = True
        if doCompress:
            startStr = outStr
            sObj = StringIO()
            gzObj = gzip.GzipFile('', 'wb', fileobj=sObj)
            gzObj.write(startStr)
            gzObj.close()
            sObj.seek(0)
            outStr = sObj.read()
            sObj.close()
            outStr = 'GzipStartTrialData' + outStr + 'GzipEndTrialData'
            compressStr = '; before comp: %sb' % len(startStr)
        else:
            compressStr = ''

        wB = sp.write(outStr) # works in background 

        nCodesSent = nCodesSent+1
        
        print('Tr %5d; serial %5db: current vals and event stream%s'
              % (nCodesSent,wB, compressStr))


    
        sys.stdout.flush()
    

def cbParseVariables(evts):
    global client

    #sys.stdout.write('cbParseVariables: found an event: code = %s, data = %r\n' %
    #                 (evts.code, evts.data))
    #sys.stdout.flush()
    if state.codecCopy is None:
        print 'ERROR - codec copy not completed'
        return

    state.clientAccessLock.acquire()


    #name = client.codec[evts.code]
    name = state.codecCopy[evts.code]
    state.clientAccessLock.release()

    # special case some names, without normal parsing
    if name == '#stimDisplayUpdate':
        for tDict in evts.data:
            if (tDict is not None
                and tDict['name'] == 'driftStimulus'):
                state.appendToEventDict(evts.code, 'driftStimulus-draw', evts.time)
                assert(tDict['action'] == 'draw')
    else:
        # default handler, state locks inside this
        state.saveVariableValue(name, evts)

        # add to event stream
        state.appendToEventDict(evts.code, evts.data, evts.time)

        # special case some AFTER normal parsing
        if name == 'experimentXmlTrialId':
            if not state.requiredExptId == evts.data:
                raise RuntimeError, 'Incorrect ExperimentId: expected %d, got %d' % (state.requiredExptId, evts.data)
            state.experimentXmlTrialId = evts.data

        

            
def sendObjectToLaser(tO, portNum):
    """Returns: number of bytes sent"""
    # Create a socket (SOCK_STREAM means a TCP socket)
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(0.1) # no more than 100ms slop
    
    # compute length etc
    data = p.dumps(tO)
    pLen = struct.pack('>L', len(data))

    # Connect to server and send data
    sock.connect((HOST, portNum))
    sock.send(pLen)  # first 4-byte length 
    sock.send(data)  # then pickled data

    # Receive data from the server and shut down
    #received = sock.recv(1024)
    sock.close()

    return (len(data)+4)

def sendLaserParams(sendDict):
    inD = state.variableCurrValues

    if inD['doBlock2SecondLaser']==1 and inD['doBlock2']==1 and inD['tBlock2TrialNumber'] == 1:
        # block2 on and this is a block2 trial, use the second laser
        tPortNum = PORT2  # global
    else:
        tPortNum = PORT  # global, default value

    tLen = sendObjectToLaser(sendDict, tPortNum)

    d = sendDict
    if d['doLinearRamp']:
        shapeStr = ('ramp %d + %d const ms long'
                    % (d['rampLengthMs'], d['rampExtraConstantLengthMs']) )
        udStr = d['rampRampDownMs']
    elif d['doPulseTrain']:
        if type(d['pulsePeriodMs']) == str:
            periodStr = 'period %s' % d['pulsePeriodMs']
        else:
             periodStr = 'period %d ms' % d['pulsePeriodMs']   
        shapeStr = ('train len %d, pulse %d, %s'
                    % (d['trainLengthMs'], d['pulseLengthMs'], periodStr))
        udStr = d['trainRampUpMs']  # same as down for now
    else:
        raise RuntimeError, 'Neither ramp nor train set'

    
    print ("Sent %db: power %7g (off %6g)mW, up/down %dms, %s"
           % (tLen, 
              chop(sendDict['pulsePowerMw'],2),
              chop(sendDict['offPowerMw'],2),
              chop(udStr,2), # same as others now
              shapeStr))
    sys.stdout.flush()


def handle_event(event):
    sys.stdout.write('Got an event: code = %s, data = %r\n' %
                     (event.code, event.data))
    sys.stdout.flush()


################################################################
## setup fns

def initializeBridge(argv, state):
    """function called by bridge script in main() to setup link to MWorks
    argv is sys.argv
    Returns: clientObj (from IPCClientConduit())
    """

    if len(argv) > 1:
        # Client-side conduit: resource name is a script argument
        conduit_resource_name = argv[1]
        #conduit_resource_name = 'python_bridge_plugin_conduit'
    else:
        # Server-side conduit: resource name is set in the experiment
        conduit_resource_name = 'server_conduit'
        print conduit_resource_name

    state.clientAccessLock.acquire()
    client = IPCClientConduit(conduit_resource_name)
    client.initialize()
    state.clientAccessLock.release()

    print "Done setting up conduit: resource %s" % conduit_resource_name
    sys.stdout.flush()

    return client


def initializeSerial():
    """function called by bridge script in main() to setup link to Blackrock through serial
    Returns: serialObj (from serial.Serial)
    """
    global sp
    
    if doTransmitCodesToBlackrock:
        portName = '/dev/cu.usbserial-A9007DvN'
        sp = serial.Serial(port=portName, 
                           baudrate=115200,
                           bytesize=serial.EIGHTBITS,
                           parity=serial.PARITY_NONE,
                           stopbits=serial.STOPBITS_ONE,
                           timeout=5.0,
                           writeTimeout=5.0,
                           xonxoff=False,
                           rtscts=False,
                           dsrdtr=False)
        print "Done setting up serial output"
        sys.stdout.flush()
    else:
        sp = None

    return sp


def initializeState(exptIdStr):
    """Set up globals"""
    global state
    state = stateC(exptIdStr)
    return state




def initializeStandardVariableCallbacks(client, state):
    """function called by bridge script in main() to setup link to Blackrock through serial
    client: clientObj (from IPCClientConduit)
    state: state obj
    Returns: None
    """
    try:
        state.clientAccessLock.acquire()

        client.register_callback_for_name('sendLaserParams', cbSendLaserParams)
        if doTransmitCodesToBlackrock:
            client.register_callback_for_name('sendSerialParams', cbSendSerialParams)
            client.register_callback_for_name('strobedDigitalWord', cbSendDigitalWordSerial)
        
        skipList = [ 'sendLaserParams',
                     'sendSerialParams',
                     'debuggerStep',
                     'debuggerRunning',
                     'strobedDigitalWord',
                     'FIO0', 'FIO1']
        keepList = [ '#stimDisplayUpdate' ]

        print "Registering callbacks, please wait...",
        sys.stdout.flush()
        
        tN = 0

        #todo: make this a method
        state.codecCopy = copy.deepcopy(client.codec)

        for name in client.codec.values():
            #print name; sys.stdout.flush()
            if ( (not (name in keepList)) 
                 and (name[0] == '#' or name in skipList) ):
                # internal event
                #print 'skipping %s' % name
                continue
            else:
                client.register_callback_for_name(name, cbParseVariables)
                tN = tN+1

        state.nCodes = tN

    finally:
        state.clientAccessLock.release()

    print (" %d codes; exptId %s (%d)" 
           % (tN, cs.reverseExperimentXmlTrialIds[state.requiredExptId], state.requiredExptId))
    print "** Ready to start **"
    sys.stdout.flush()

    #import rpdb2; rpdb2.start_embedded_debugger('password')



def finalizeBridge(client, state):
    """MH 130102: not sure when this ever gets run"""

    print "Finalizing conduit...",
    sys.stdout.flush()
    if not state.clientAccessLock.locked():
        state.clientAccessLock.acquire()
    client.finalize()
    state.clientAccessLock.release()
    print "done"
    sys.stdout.flush()

