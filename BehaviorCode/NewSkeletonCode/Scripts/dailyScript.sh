#!/bin/sh
#  dailyScript.sh
#  
# Depends on RSA handshake between HullGlickN computers.
#  Created by Andrew McKinney on 11/8/2014.
#
scp -r hullglick@hullglick2.local:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick2.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
scp -r hullglick@hullglick3.local:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick3.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
scp -r hullglick@hullglick4.local:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick4.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks

mount -t smbfs //andrew:dukegrad@crash.dhe.duke.edu/andrew/behavior ~/Desktop/Crash_Server
scp -r --exclude=.DS_Store ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
scp -r --exclude=.DS_Store ~/Documents/MWorks/Data ~/Desktop/Crash_server/

#/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "dailyScript([212 213 214 305 306]);"
#/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "quit"
#scp -r ~/Documents/MWorks/DailyPlots ~/Desktop/Crash_server/

umount ~/Desktop/Crash_Server

exit
