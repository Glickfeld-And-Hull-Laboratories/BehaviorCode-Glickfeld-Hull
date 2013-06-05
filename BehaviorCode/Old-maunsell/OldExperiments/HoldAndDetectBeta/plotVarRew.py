# variable reward plot if using a timing strategy

import matplotlib.pyplot as plt
import numpy as np

class v:
    pass   #make a struct

v.randReqHoldMaxMs = 900
v.fixedReqHoldTimeMs = 100
v.tooFastTimeMs = 1
v.fixedEarlyTimeMs = 150
v.earlyRewardUs = 3000
v.minRewardUs = 30000
v.maxRewardUs = 45000
v.minAgainRewardUs = 3000

v.varRewActualMinMs = 100
v.varRewActualMaxMs = 250
v.varRewActualMinAgainMs = 700
v.reactTimeMs = 3000

#figH = plt.figure()


nRuns = 2000


def doSim(nRuns, doVis, v):

    nPts = 4000
    allVR = np.zeros((nPts, nRuns))

    for iR in range(0,nRuns):  
        # simulate nRuns trials
        tVarRew = np.zeros((nPts,))
    
        if doVis:
            tHoldTime = (v.randReqHoldMaxMs*np.random.random_sample(1) 
                         + v.fixedReqHoldTimeMs)
            t = np.asarray(np.round(tHoldTime), dtype='int64')
        else:
            error()



        tVarRew[0:v.fixedEarlyTimeMs] = v.earlyRewardUs
        tVarRew[v.fixedEarlyTimeMs:t] = v.earlyRewardUs
        tVarRew[t:t+v.varRewActualMinMs] = v.minRewardUs
        tVarRew[t+v.varRewActualMinMs:t+v.varRewActualMaxMs] \
            = np.linspace(v.minRewardUs, v.maxRewardUs, 
                          v.varRewActualMaxMs-v.varRewActualMinMs)
        tVarRew[t+v.varRewActualMaxMs:t+v.varRewActualMinAgainMs] \
            = np.linspace(v.maxRewardUs, v.minAgainRewardUs, 
                          v.varRewActualMinAgainMs-v.varRewActualMaxMs)
        tVarRew[t+v.varRewActualMinAgainMs:t+v.reactTimeMs] = v.minAgainRewardUs
        tVarRew[t:t+v.tooFastTimeMs] = v.earlyRewardUs
 
       #else not written code yet
        allVR[:,iR] = tVarRew

    return allVR

allVR = doSim(nRuns, 1, v)
meanVR = np.mean(allVR, axis=1)

plt.plot(meanVR)

t1 = np.round(v.randReqHoldMaxMs/2 + v.fixedReqHoldTimeMs)
plt.plot((0,t1,t1, t1+v.reactTimeMs), (0, 0, v.maxRewardUs, v.maxRewardUs))

plt.gca().set_xlim(0,2500)




plt.title('Max %d ms' % np.argmax(meanVR))

# wait when run from command line
#canvas = plt.gcf().canvas
#canvas.start_event_loop(timeout=0)

