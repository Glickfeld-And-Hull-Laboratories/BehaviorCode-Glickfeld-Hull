import signal
import sys
import socket
import sys
import cPickle as p
import struct
import time
import numpy as np
import math
import serial
import json





sp2 = serial.Serial(port='/dev/cu.usbserial-A9007DvM',
                    baudrate=115200,
                    bytesize=serial.EIGHTBITS,
                    parity=serial.PARITY_NONE,
                    stopbits=serial.STOPBITS_ONE,
                    timeout=5.0,
                    writeTimeout=5.0,
                    xonxoff=False,
                    rtscts=False,
                    dsrdtr=False)

testDict = { 'a': 1, 'b': 2 }

tStr = json.dumps(testDict);

print tStr;

for iR in range(100):
    sp.write(tStr);
    sp.write('\n');
    time.sleep(1);
    sp2.write(tStr);
    sp.write('\n');
    time.sleep(1)
    
    nToRead = sp2.inWaiting()
    outB = sp2.read(nToRead)
    print '%d: %s' % (nToRead, outB)

    nToRead = sp.inWaiting()
    outB = sp.read(nToRead)
    print '%d: %s' % (nToRead, outB)
    
    time.sleep(1)
    


                   
