#!/bin/sh
#  dailyScript.sh
#  
# Depends on RSA handshake between HullGlickN computers.
#  Created by Andrew McKinney on 11/12/13.
#

#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick2.local:~/Documents/MWorks/Data ~/Documents/MWorks 
#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick2.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick3.local:~/Documents/MWorks/Data ~/Documents/MWorks
#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick3.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks
#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick4.local:~/Documents/MWorks/Data ~/Documents/MWorks
#sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick4.local:~/Documents/MWorks/BehavOutputPdfs ~/Documents/MWorks

mount -t smbfs //andrew:dukegrad@crash.dhe.duke.edu/andrew/behavior ~/Desktop/Crash_Server
for num in 1 2 3 4
do
    sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick$num.local:~/Documents/MWorks/Data ~/Desktop/Crash_server/ 
    sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick$num.local:~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
done

#sshpass -p hullglick rsync --ignore-existing -rv ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
#sshpass -p hullglick rsync --ignore-existing -rv ~/Documents/MWorks/Data ~/Desktop/Crash_server/

#/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "dailyScript([212 213 214 305 306]);"
#/Applications/MATLAB_R2012b.app/bin/matlab -nosplash -nodesktop -r "quit"
#rsync -rv ~/Documents/MWorks/DailyPlots ~/Desktop/Crash_server/

umount ~/Desktop/Crash_Server

exit