import signal
import sys
import socket
import sys
import cPickle as p
import struct
import time
import numpy as np
import math

from mworks.conduit import IPCClientConduit

sentOnce = False

def handle_event(event):
    sys.stdout.write('Got an event: code = %s, data = %r\n' %
                     (event.code, event.data))

    if not sentOnce:
        print "calling client.send_etc"
        client.send_integer(100, np.long(41))
        # also - using the code from the codec (in this case 81) does not work
        # neither does send_data or send_float
        #client.send_data(100,41)
        sys.stdout.flush()
        sentOnce = True


def main():
    global client
    
    conduit_resource_name = sys.argv[1]
    client = IPCClientConduit(conduit_resource_name)
    client.initialize()
    print "Done setting up"; sys.stdout.flush()

    try:
        client.register_callback_for_name('ignore', handle_event)
        client.register_local_event_code(100, 'ignore')

        client.send_float(100, 40.0)  # no effect

        print "Registered"; sys.stdout.flush()
        
        signal.pause()  # Other threads are working; this one can sleep
    finally:
        client.finalize()

    print "End"
    sys.stdout.flush()


if __name__ == '__main__':
    main()
