import numpy as np

def chopUpFloatIntoSigArray(inFloat, nSigDigits):
    """Returns: np.array(dtype=uint8) of length nSigDigits+1.  The last byte is the exponent.
    If number is negative first entry will be negative (two's complement)
    This algorithm is slow as heck but the code should be clear """
    outArr = np.zeros(nSigDigits+1, dtype='uint8')

    # deal with negative numbers up front
    if inFloat < 0:
        firstSign = -1
        inFloat = inFloat * -1
    else:
        firstSign = 1

    workN = inFloat
    for iB in np.arange(nSigDigits):
        tD = 10**np.floor(np.log10(workN))
        tSigD = (np.uint8)(np.trunc(workN/tD))
        workN = workN - tSigD*tD
        if iB == 0 and firstSign == -1:
            outArr[iB] = -tSigD
        else:
            outArr[iB] = tSigD
            
    outArr[nSigDigits] = np.floor(np.log10(inFloat))
    
    return outArr

b=chopUpFloatIntoSigArray(122342, 3)
