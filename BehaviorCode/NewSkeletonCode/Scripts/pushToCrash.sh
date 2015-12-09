mount -t smbfs //andrew:dukegrad@crash.dhe.duke.edu/andrew/behavior ~/Desktop/Crash_Server
rsync --ignore-existing -rv ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
rsync --ignore-existing -rv ~/Documents/MWorks/Data ~/Desktop/Crash_server/
umount ~/Desktop/Crash_Server

exit
