#!/bin/bash
# push-files-one-rig.sh MaunsellMouseN.local [doOneShot]

### switch these for debug mode

#RSYNC_OPTS="-av --dry-run"
RSYNC_OPTS="-av"

if [[ "$2" == "1" ]]; then
    DO_ONE_SHOT=1
else
    DO_ONE_SHOT=0
fi

####

outcomputer=$1
TARGUSER=holdanddetect

shopt -s extglob
set -e  # exit on error, esp on the scp below

# Make sure we have a Snow Leopard install on the target disk
echo -n Looking for remote: Snow Leopard:
scp -rpq -B "$outcomputer:/Library/Application\ Support/MWorks/Configuration/setup_variables.xml" /tmp

echo -n found, user: $TARGUSER:
scp -rpq -B "$outcomputer:/Users/$TARGUSER/.ssh/known_hosts" /tmp 
echo found

# Applications
echo -n . MWorks applications...
rsync $RSYNC_OPTS --delete /Applications/?(MWClient.app|MWServer.app|MWEditor.app) root@$outcomputer:/Applications >> /tmp/push-log.txt 
echo done

# /usr/local tree  
# dangerous!  do not do this - mysql/data server problems, as well as bin/
#rsync $RSYNC_OPTS --delete --dry-run --exclude 'mysql*' /usr/local/* root@$outcomputer:/usr/local

# more limited
#usrLocalOpts="--dry-run"
usrLocalOpts=""
echo -n . Parts of /usr/local...
rsync $RSYNC_OPTS --delete $usrLocalOpts --exclude 'mysql*' /usr/local/lib root@$outcomputer:/usr/local >> /tmp/push-log.txt 
rsync $RSYNC_OPTS --delete $usrLocalOpts --exclude 'mysql*' /usr/local/include root@$outcomputer:/usr/local >> /tmp/push-log.txt 
rsync $RSYNC_OPTS --delete $usrLocalOpts --exclude 'mysql*' /usr/local/Library root@$outcomputer:/usr/local >> /tmp/push-log.txt 
echo done

# Frameworks
echo -n . /Library/Frameworks...
rsync $RSYNC_OPTS --delete /Library/Frameworks/?(MWorksCore|MWorksCocoa|Sparkle).framework "root@$outcomputer:/Library/Frameworks" >> /tmp/push-log.txt 
echo done

# App support
echo -n . /Library/Application Support ...
rsync $RSYNC_OPTS --delete --exclude MATLAB-old "/Library/Application Support/MWorks" "root@$outcomputer:/Library/Application\ Support" >> /tmp/push-log.txt 
echo done

# expt xml
echo -n . ExperimentXML-git/ and Repositories/ ...
rsync $RSYNC_OPTS --delete --exclude generated-images /Users/$TARGUSER/ExperimentXML-git root@$outcomputer:/Users/$TARGUSER >> /tmp/push-log.txt 
#rsync $RSYNC_OPTS --delete /Users/$TARGUSER/Repositories/MWorksMatlabToolbox root@$outcomputer:/Users/$TARGUSER/Repositories >> /tmp/push-log.txt 
rsync $RSYNC_OPTS --delete --exclude MWorks /Users/$TARGUSER/Repositories root@$outcomputer:/Users/$TARGUSER >> /tmp/push-log.txt 
#rsync $RSYNC_OPTS --delete /Users/$TARGUSER/ExperimentXML-Intervals root@$outcomputer:/Users/$TARGUSER >> /tmp/push-log.txt 
echo done

# various configs
#rsync $RSYNC_OPTS /Users/$TARGUSER/.ssh/authorized_keys root@$outcomputer:/Users/$TARGUSER/.ssh/authorized_keys
echo -n . dotfiles ...
rsync $RSYNC_OPTS /Users/$TARGUSER/.bash_profile root@$outcomputer:/Users/$TARGUSER/ >> /tmp/push-log.txt 
echo done


# one time shots

#rsync $RSYNC_OPTS --delete /Applications/MATLAB_R2008b.app root@$outcomputer:/Applications
#rsync $RSYNC_OPTS --delete /Applications/MATLAB root@$outcomputer:/Applications
#rsync $RSYNC_OPTS --delete /Library/Frameworks/Python.framework root@$outcomputer:/Library/Frameworks
#rsync $RSYNC_OPTS --delete /Users/$TARGUSER/MWorksSoftware root@$outcomputer:/Users/$TARGUSER
#rsync $RSYNC_OPTS --delete /Applications/Aquamacs* root@$outcomputer:/Applications
#rsync $RSYNC_OPTS --delete /Applications/Tower.app root@$outcomputer:/Applications
if [ $DO_ONE_SHOT -eq 1 ]; then
    echo
    echo "******* Doing one-shot syncs **********"
    echo
    rsync $RSYNC_OPTS --delete --exclude \*.lic /Applications/MATLAB_R2011a.app root@$outcomputer:/Applications
    rsync $RSYNC_OPTS --delete /Applications/Aquamacs* root@$outcomputer:/Applications
    rsync $RSYNC_OPTS /Users/$TARGUSER/.git* root@$outcomputer:/Users/$TARGUSER/
    rsync $RSYNC_OPTS --delete /Users/$TARGUSER/Repositories root@$outcomputer:/Users/$TARGUSER
    rsync $RSYNC_OPTS --delete /Applications/Tower.app root@$outcomputer:/Applications
fi


