import math
import numpy
import random


imageList = [] #create empty list for images desired
setsWanted = getvar('setList') #get list of names of image set

#add to this when new image set is desired
if 'flowers' in setsWanted: 
	imageList.extend([1, 2, 3, 4, 5, 6, 7])
if 'fruits' in setsWanted:
	imageList.extend([8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
if 'landscape' in setsWanted:
	imageList.extend([24, 25, 26, 27, 28, 29, 30])



setvar('stimList', imageList) #set list of images desired in mworks
