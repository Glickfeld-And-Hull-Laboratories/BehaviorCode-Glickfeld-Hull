import sys
import signal
import ExptBase
import copy
import ExptConstants as cs
import numpy as np

global varVects
varVects = {}

# numpy useful utils
def futureIn1d(a,b):
    """delete this when we have numpy 1.4"""
    return np.array([item in b for item in a])

def cbSync(event):
    global state
    global client
    
    
    if event.data == 1:  
        #print 'start'; sys.stdout.flush()
        state.resetThisTrialData()

    elif event.data == 0:
        #print 'end'
        if state.variableThisTrialValues.has_key('subjectNum'):
            #print 'Constant dump found'
            #sys.stdout.flush()
            # now do nothing with constants, they get updated automatically into variableCurrValues
            pass
        else:
            #print 'found real end of trial'
            #sys.stdout.flush()
            cbEndOfTrial()

    else:
        raise RuntimeError

def listUpdateIndex(tList,index,value):
    """ updates in place """
    if (index+1) > len(tList):
        tList.extend([None]*(index-len(tList)+1))
    tList[index] = value
    
def cbEndOfTrial():
    global state
    global client
    global varVects
    debug = True

    ## save vector data: one value each trial ################

    varD = copy.copy(state.variableThisTrialValues)
    constD = state.getVariableValues()
    nTrs = varD['tNTrialsCompleted']
    
    # init on first trial
    vectList = ['outcomeList', 'tLeftTrial']
    for vN in vectList:
        if not varVects.has_key(vN):
            varVects[vN] = [] # init to empty list


    # do trial outcomes
    assert (not varD.has_key('subjectNum')), 'should never be called for constant dump: bug'
    npTrialOutcomes = np.array(cs.trialOutcomes)
    foundIx = futureIn1d(npTrialOutcomes, varD.keys())
    if sum(foundIx) != 1:
        print foundIx
        print varD.keys()
        raise RuntimeError, 'no outcome found for this trial: bug'
    tOutcomeStr = npTrialOutcomes[foundIx][0]
    listUpdateIndex(varVects['outcomeList'], nTrs-1, tOutcomeStr)

    # update single element lists
    listUpdateIndex(varVects['tLeftTrial'], nTrs-1, varD['tLeftTrial'])

    ## do the bias computation ################
    
    # consts
    biasNPts = constD['rewardStaircaseWindowWidth']
    biasMinTrs = constD['rewardStaircaseWindowWidth']
    nTrToUse = np.min((len(varVects['outcomeList']),biasNPts))

    # get data vectors, and compute response/left vectors
    varArr = {}
    for n in varVects.keys():
        varArr[n] = np.array(varVects[n])

    # vectors of correct and incorrect trials as used by the staircase
    stairCorrV = np.array(futureIn1d(varArr['outcomeList'], ['success']), dtype='bool')
    stairNotCorrV = np.array(futureIn1d(varArr['outcomeList'], ['incorrect']), dtype='bool')

    # create a mask vector giving all the trials to consider (in the window)
    considerV = (stairCorrV == 1) | (stairNotCorrV == 1)
    accumNFromEnd = np.nonzero(np.cumsum(considerV[::-1])==nTrToUse)[0]  # reverse, sum bool to match count
    if len(accumNFromEnd) > 0:
        firstToConsiderN = len(considerV) - accumNFromEnd[0] -1  # , subtr from len to get index from start
    if len(accumNFromEnd) <= 0 or firstToConsiderN < 0:
        firstToConsiderN = 0
    trToConsiderIx = np.zeros(nTrs, dtype='bool')
    trToConsiderIx[firstToConsiderN:] = True
    nToConsider = np.sum(trToConsiderIx)

    if debug:
        print '%d %d %d %d' % (nTrToUse, len(accumNFromEnd), firstToConsiderN, len(trToConsiderIx))

    
    if constD['rewardStaircaseDoEqualizeCorrects']:
        leftCorrN = np.sum(stairCorrV[ (varArr['tLeftTrial']==1) & (trToConsiderIx == 1) ] == 1)
        leftNotCorrN = np.sum(stairNotCorrV[ (varArr['tLeftTrial']==1) & (trToConsiderIx == 1) ] == 1)
        rightCorrN = np.sum(stairCorrV[ (varArr['tLeftTrial']==0) & (trToConsiderIx == 1) ] == 1)
        rightNotCorrN = np.sum(stairNotCorrV[ (varArr['tLeftTrial']==0) & (trToConsiderIx == 1) ] == 1)
        leftCorrFract = leftCorrN / float(leftCorrN + leftNotCorrN)
        rightCorrFract = rightCorrN / float(rightCorrN + rightNotCorrN)
        computedBias = (leftCorrFract - rightCorrFract) / 2 + 0.5 # map range from -1..+1 to 0..1
        if debug:
            print ('left %d %d %3.2f; right %d %d %3.2f; bias %3.2f' 
                   % (leftCorrN, leftNotCorrN, leftCorrFract, 
                      rightCorrN, rightNotCorrN, rightCorrFract, 
                      computedBias))
        
    else:  
        # equalize responses each side
        respDirLeft = np.array([None]*nTrs, dtype='float64')
        respDirLeft[np.logical_and(varArr['tLeftTrial']==1, varArr['outcomeList'] == 'success')] = 1
        respDirLeft[np.logical_and(varArr['tLeftTrial']==0, varArr['outcomeList'] == 'incorrect')] = 1
        respDirLeft[np.logical_and(varArr['tLeftTrial']==0, varArr['outcomeList'] == 'success')] = 0
        respDirLeft[np.logical_and(varArr['tLeftTrial']==1, varArr['outcomeList'] == 'incorrect')] = 0

        # compute bias
        respDirOnlyResp = respDirLeft[~np.isnan(respDirLeft)]
        nTrToUse = np.min((len(varVects['outcomeList']),biasNPts))
        computedBias = np.sum(respDirOnlyResp[-nTrToUse:])/nTrToUse

    # correct bias
    if nTrs > biasMinTrs or np.isnan(computedBias):
        if np.isnan(computedBias): 
            assert (firstToConsiderN == 0)  # else computed bias is nan after many trials are run: bug
        tWinLeftBias = computedBias
    else:
        tWinLeftBias = 0.5  # will be a constant value if before min trials

    # last check
    if np.isnan(tWinLeftBias):
        print '** BUG **: tWinLeftBias is nan, check and fix'
        # continue
        tWinLeftBias = constD['rewardStaircaseStartingLeftBias']
        #tWinLeftBias = 0.5  
        
    # write bias to MWorks, on every trial, if DoSlidingWin is true
    if constD['rewardStaircaseDoSlidingWin'] == 1:
        client.send_data(client.reverse_codec['tLeftBias'], tWinLeftBias)
        wroteStr = 'sent to MWorks'
    else:
        wroteStr = 'not sent to MWorks'

    print 'sliding win %d resp trials (%d total); tRunningLeftBias: %4.2f, %s' % (nTrToUse, nToConsider, tWinLeftBias, wroteStr)
    sys.stdout.flush()





    

    
    
def main():
    global sp
    global client
    global state  # ini
    exptIdStr = 'DualPress'


    try:
        state = ExptBase.initializeState(exptIdStr)
        client = ExptBase.initializeBridge(sys.argv, state)
        sp = ExptBase.initializeSerial()
        ExptBase.initializeStandardVariableCallbacks(client, state)

        # register a end-of-trial-callback
        client.register_callback_for_name('sync', cbSync)

        signal.pause()  # Other threads are working; this one can sleep

    finally:
        ExptBase.finalizeBridge(client, state)

if __name__ == '__main__':
    main()


# to convert to HADC8:
# need a tNTrialsCompleted variable
