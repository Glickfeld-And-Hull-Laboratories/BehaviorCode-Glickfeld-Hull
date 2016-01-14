import sys
import signal
import ExptBase


def main():
    global sp
    global client
    global state  # ini
    exptIdStr = 'HoldAndDetectConstant8'

    try:
        state = ExptBase.initializeState(exptIdStr)
        client = ExptBase.initializeBridge(sys.argv, state)
        sp = ExptBase.initializeSerial()
        ExptBase.initializeStandardVariableCallbacks(client, state)
        print('dump')

        
        signal.pause()  # Other threads are working; this one can sleep

    finally:
        pass
        #ExptBase.finalizeBridge(client, state)

if __name__ == '__main__':
    main()
