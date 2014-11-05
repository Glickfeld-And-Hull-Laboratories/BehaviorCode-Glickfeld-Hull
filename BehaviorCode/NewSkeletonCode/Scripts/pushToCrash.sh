mount -t smbfs //andrew:dukegrad@crash.dhe.duke.edu/andrew/behavior ~/Desktop/Crash_Server
scp -r ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Crash_server/
scp -r ~/Documents/MWorks/Data ~/Desktop/Crash_server/
umount ~/Desktop/Crash_Server

exit
