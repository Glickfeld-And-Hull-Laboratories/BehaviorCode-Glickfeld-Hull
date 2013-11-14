#!/bin/sh
#  dailyScript.sh
#  
# Depends on RSA handshake between HullGlickN computers.
#  Created by Andrew McKinney on 11/12/13.
#
scp -r hullglick@hullglick2.dhe.duke.edu:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick2.dhe.duke.edu:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
scp -r hullglick@hullglick3.dhe.duke.edu:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick3.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
scp -r hullglick@hullglick4.local:~/Documents/MWorks/Data ~/Documents/MWorks
scp -r hullglick@hullglick4.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks

mount -t smbfs //andrew:dukegrad@crash.dhe.duke.edu/andrew/behavior ~/Desktop/Crash_Server
scp -r ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
scp -r ~/Documents/MWorks/Data ~/Desktop/Crash_server/

/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "dailyScript([1 5 6 202 205 206 209]);"
/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "quit"
scp -r ~/Documents/MWorks/DailyPlots ~/Desktop/Crash_server/

umount ~/Desktop/Crash_Server

exit
