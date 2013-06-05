import socket
import cPickle as p
import struct
import sys
import time
import random

HOST = 'mark-laser-PC.local'
PORT = 9990

laserOutDict = {'pulsePowerMw':10 , 'offPowerMw':5, 'pulseLengthMs':100}

def sendObject(sO):
    """returns:
    number of bytes sent
    """
    
    # Create a socket (SOCK_STREAM means a TCP socket)
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(1) # no more than 100ms slop
    
    # compute length etc
    data = p.dumps(sO)
    pLen = struct.pack('>L', len(data))

    # Connect to server and send data
    sock.connect((HOST, PORT))
    sock.send(pLen)  # first 4-byte length 
    sock.send(data)  # then pickled data

    # Receive data from the server and shut down
    #received = sock.recv(1024)
    sock.close()

    return(len(data)+4)

def sendLaserParams(sendDict):

    dLen = sendObject(sendDict)

    print ("Sent: %d bytes; power %6g-%6g, length %dms"
           % (dLen, sendDict['offPowerMw'], sendDict['pulsePowerMw'], sendDict['pulseLengthMs']))
    sys.stdout.flush()


while True:
    # print "Trying to send..."
    # try:
    #     sendLaserParams(laserOutDict)
    # except:
    #     pass

    sendLaserParams(laserOutDict)
    time.sleep(0.2)
    
    dLen = sendObject('abortOutput')
    print ("Sent: %d bytes" % dLen)
    time.sleep(1)
            
