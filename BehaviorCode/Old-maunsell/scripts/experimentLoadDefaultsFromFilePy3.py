#!/usr/bin/python3
#  requires python2.5; tested only on Snow Leopard, will probably not work on Lion or later.  
#    (for upgrade: check the PyObjC calls
#
#  first argument: file to read, defaults to ~/Desktop/org.histed.MWClientSavedVars.plist
#
#  histed 121022

from __future__ import with_statement

import objc
from Foundation import (
	NSMutableDictionary, 
	NSUserDefaults, 
	NSDictionary, 
	NSCFArray, 
)
from copy import copy
import os
import re
import sys

domainName = "org.mworks-project.MWClient"
homedir = os.getenv('HOME')

################

def subStr(inStr):
    return re.sub('^\\$HOME', homedir, inStr)

def replaceStrWithUserdir(inObj):
    if type(inObj) == str or type(inObj) == objc.pyobjc_unicode:
        return subStr(inObj)
    elif isinstance(inObj, NSCFArray):
        for i in range(len(inObj)):
            # do this recursively
            inObj[i] = replaceStrWithUserdir(inObj[i])
        return inObj
    else:
        print (type(inObj))
        #import pdb; pdb.set_trace()
        raise Exception('Error: Type unknown')
    return 

################

readFile = sys.argv[1]
if readFile is None:
  readFile = os.path.join(os.getenv('HOME'), 'Desktop', 'org.histed.MWClientSavedVars.plist')

# get defaults from file
fileDefs = NSDictionary.dictionaryWithContentsOfFile_(readFile)

# get user defaults objects
standardUserDefaults = NSUserDefaults.standardUserDefaults()
clientDefs = standardUserDefaults.persistentDomainForName_(domainName).mutableCopy()

keyList = fileDefs.keys()
for k in keyList:
  # replace $HOME strs
  #fileDefs[k] = replaceStrWithUserdir(fileDefs[k])
  clientDefs.setObject_forKey_(fileDefs[k], k)

  
standardUserDefaults.setPersistentDomain_forName_(clientDefs, domainName)

standardUserDefaults.synchronize()



