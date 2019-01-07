mount -t smbfs //jarvis_lab@duhs-user-nc1.dhe.duke.edu/dusom_glickfeldlab/All_Staff/Behavior ~/Desktop/Isilon_Server
rsync --ignore-existing -rv ~/Documents/MWorks/BehavOutputPdfs ~/Desktop/Isilon_server/
rsync --ignore-existing -rv ~/Documents/MWorks/Data ~/Desktop/Isilon_server/

umount /Volumes/Behavior
exit
