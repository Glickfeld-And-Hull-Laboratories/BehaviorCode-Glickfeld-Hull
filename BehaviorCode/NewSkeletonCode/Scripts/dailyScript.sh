#!/bin/sh
#  dailyScript.sh
#  
# Depends on RSA handshake between HullGlickN computers.
#  Created by Andrew McKinney on 11/12/13.
# Updated by LG on 11/30/18 to connect to Isilon instead (using generic jarvis_lab account)

mount -t smbfs //jarvis_lab@duhs-user-nc1.dhe.duke.edu/dusom_glickfeldlab/All_Staff/Behavior ~/Desktop/Isilon_Server
for num in 1 2 3 4
do
    sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick$num.local:~/Documents/MWorks/Data ~/Desktop/Isilon_server/
    sshpass -p hullglick rsync --ignore-existing -rv hullglick@hullglick$num.local:~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Isilon_server/
done


umount ~/Desktop/Isilon_Server

exit
