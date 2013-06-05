import u6
import sys, os, time
from datetime import datetime

d = u6.U6()

timesOpenMs = [75,   15,  20,  25,  30,  35,  40,  50,  75,    5]
nRepsMs     = [500, 500, 500, 500, 500, 500, 500, 500, 500, 1000]
closedBetweenMs = 50;

u6IONumber = 0 # FIO0

try:

    d.getFeedback(u6.BitDirWrite(u6IONumber, 1)) # output

    print "waterCalib.py:  Waiting for %d ms between pulses" % closedBetweenMs

    for (timeOpen,nReps) in zip(timesOpenMs, nRepsMs):
        print ("Doing %d ms data point, %d reps; set up catcher and press return"
               % (timeOpen, nReps))
        raw_input()

        timeStartT = time.time()
        for iR in range(0,nReps):
            startT = time.time()

            # open
            d.getFeedback(u6.BitStateWrite(u6IONumber, 1)) # high

            time.sleep( ((timeOpen-5)/1000) ) # sleep for all but 5ms
        
            while (time.time()-startT) * 1000 <= timeOpen:
                # spin and wait for timer to expire
                pass

            d.getFeedback(u6.BitStateWrite(u6IONumber, 0)) # low

            time.sleep(closedBetweenMs/1000.0)

            if iR % 50 == 0:
                sys.stdout.write("  %d" % iR)
                sys.stdout.flush()

        sys.stdout.write("\n")
    
        print "Done, did %d reps in %6.3fs" % (nReps, time.time()-timeStartT)

finally:

    print "** Stopped, cleaning up"
    d.getFeedback(u6.BitStateWrite(u6IONumber, 0)) # set low
    d.close() # so next start works ok
