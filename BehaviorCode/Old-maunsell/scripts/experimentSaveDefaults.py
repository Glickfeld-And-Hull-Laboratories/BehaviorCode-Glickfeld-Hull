#!/usr/bin/python2.5
#
#  requires python2.5; tested only on Snow Leopard, will probably not work on Lion or later.  
#    (for upgrade: check the PyObjC calls
#
#  writes to a file on the Desktop by default
#  searches for home directory in each key by string matching and replaces w/ '$HOME'
#
#  histed 121022
#  edited mckinney 140115

from __future__ import with_statement

from Foundation import NSMutableDictionary, NSUserDefaults, NSCFArray, objc
import numpy as np
from copy import copy
import os
import re

domainName = "org.mworks-project.MWClient"

outFile = os.path.expanduser(os.path.join('~/Repositories/BehaviorCode-Glickfeld-Hull/BehaviorCode/NewSkeletonCode/Scripts/Plists', 'org.Behavior.flashingStim2.plist'))
keyNames = [ 
  'MATLAB client window - selected variables',
  'MATLAB client window - MATLAB .m file',
  'recentPythonScripts' ]

homedir = os.getenv('HOME')

################

def subStr(inStr):
    return re.sub('^%s'%homedir, '$HOME', inStr)

def replaceUserdirWithStr(inObj):
    if type(inObj) == str or type(inObj) == objc.pyobjc_unicode:
        return subStr(inObj)
    elif isinstance(inObj, NSCFArray):
        for i in range(len(inObj)):
            # do this recursivel#y
            inObj[i] = replaceUserdirWithStr(inObj[i])
        return inObj
    else:
        print type(inObj)
        #import pdb; pdb.set_trace()
        raise Exception('Error: Type unknown')
    return 
    

################

# get client defaults
standardUserDefaults = NSUserDefaults.standardUserDefaults()
clientDefs = standardUserDefaults.persistentDomainForName_(domainName)

# copy the fields we need
writeDict = NSMutableDictionary.dictionary()
for k in clientDefs:
  if k in keyNames:
    tVal = clientDefs[k]
    #tVal = replaceUserdirWithStr(tVal)
    writeDict[k] = tVal

# use the Cocoa call to write output as XML
success = writeDict.writeToFile_atomically_(outFile, 1)



