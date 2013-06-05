import signal
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

# set machine constants
thisHostname = socket.gethostname()
thisHostname = string.replace(thisHostname, '.local', '')
if thisHostname == "MaunsellMouse1":
    PORT = 9990
    doTransmitCodesToBlackrock = True
elif thisHostname == "MaunsellMouse2":
    PORT = 9991
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouse3":
    PORT = 9992
    doTransmitCodesToBlackrock = False
elif thisHostname == "MaunsellMouse4":
    PORT = 9993
    doTransmitCodesToBlackrock = False
else:
    raise RuntimeError, 'Found unknown hostname: %s' % thisHostname
    

HOST = "mark-laser-PC.local"
laserParamCode = 'sendLaserParams'
lastSendTime = 0
nCodesSent = 0

if doTransmitCodesToBlackrock:
    import json
    import serial

class constantsC:
    experimentXmlTrialIds = { 'HoldAndDetectConstant8' : 8,
                              'HoldAndDetectConstant5' : 10,
                              'DirTuningMapping' : 11,
                              'ChRMappingTwelve' : 12,
                              'ChRMappingTen': 13}

class stateC:
    lock = None
    nCodes = None
    nStrobedTotal = 0
    eventDict = None
    experimentXmlTrialId = None
    blankEventDict = { 'names':[], 'values':[], 'ts':[] }
    variableCurrValues = {}
    variableLastTs = {}
    clientAccessLock = None
    codecCopy = None

    def __init__(self):
        self.lock = threading.Lock()
        self.clientAccessLock = threading.Lock()

        self.eventDict = deepcopy(self.blankEventDict)

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


state = stateC()
cs = constantsC()

################

def computeLaserPower():
    """pass in variableCurrValues,
    returns a dictionary that can be passed over the wire to the laser stim machine"""
    inD = state.variableCurrValues
    
    outD = { 'pulsePowerMw': float(inD['tLaserPowerMw']),
             'offPowerMw': float(inD['laserOffPowerMw']),
             'doLinearRamp': int(inD['tLaserDoLinearRamp']),
             'rampLengthMs': float(inD['tLaserRampLengthMs']),
             'doPulseTrain': int(inD['tLaserDoPulseTrain']),
             'pulseLengthMs': float(inD['laserPulseLengthMs']),
             'pulsePeriodMs': float(inD['laserPulsePeriodMs']),
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
                laserOutDict = computeLaserPower()
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
        
        # default handler
        state.lock.acquire() # needed only during writes to data structures
        state.variableCurrValues[name] = evts.data
        state.variableLastTs[name] = evts.time
        state.lock.release()

        # add to event stream
        state.appendToEventDict(evts.code, evts.data, evts.time)

        # special case some AFTER normal parsing
        if name == 'experimentXmlTrialId':
            state.experimentXmlTrialId = evts.data

        

            
def sendObjectToLaser(tO):
    """Returns: number of bytes sent"""
    # Create a socket (SOCK_STREAM means a TCP socket)
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(0.1) # no more than 100ms slop
    
    # compute length etc
    data = p.dumps(tO)
    pLen = struct.pack('>L', len(data))

    # Connect to server and send data
    sock.connect((HOST, PORT))
    sock.send(pLen)  # first 4-byte length 
    sock.send(data)  # then pickled data

    # Receive data from the server and shut down
    #received = sock.recv(1024)
    sock.close()

    return (len(data)+4)

def sendLaserParams(sendDict):
    tLen = sendObjectToLaser(sendDict)

    d = sendDict
    if d['doLinearRamp']:
        shapeStr = ('ramp %d + %d const ms long'
                    % (d['rampLengthMs'], d['rampExtraConstantLengthMs']) )
        udStr = d['rampRampDownMs']
    elif d['doPulseTrain']:
        shapeStr = ('train len %d, pulse %d, period %d ms'
                    % (d['trainLengthMs'], d['pulseLengthMs'], d['pulsePeriodMs']))
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


def main():
    global sp
    global client

    if len(sys.argv) > 1:
        # Client-side conduit: resource name is a script argument
        conduit_resource_name = sys.argv[1]
        #conduit_resource_name = 'python_bridge_plugin_conduit'
    else:
        # Server-side conduit: resource name is set in the experiment
        conduit_resource_name = 'server_conduit'
        print conduit_resource_name

    state.clientAccessLock.acquire()
    client = IPCClientConduit(conduit_resource_name)
    client.initialize()

    print "Done setting up conduit: resource %s" % conduit_resource_name
    sys.stdout.flush()

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


    try:
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

            #if np.mod(tN,10) == 0:
            #    print ("%d" % tN),
            #   sys.stdout.flush()



        state.clientAccessLock.release()
        state.nCodes = tN

        print " %d codes" % (tN)
        print "** Ready to start **"
        sys.stdout.flush()

        #import rpdb2; rpdb2.start_embedded_debugger('password')

        signal.pause()  # Other threads are working; this one can sleep

    finally:
        print "Finalizing..."
        sys.stdout.flush()
        if not state.clientAccessLock.locked():
            state.clientAccessLock.acquire()
        client.finalize()
        state.clientAccessLock.release()

if __name__ == '__main__':
    main()
