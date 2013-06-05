#!/bin/bash

set -e

#thisDef=$(defaults read "org.mworks-project.MWClient" "MATLAB client window - selected variables")
#echo $thisDef
#ssh $1 defaults write "'org.mworks-project.MWClient'" "'MATLAB client window - selected variables'" \'"$thisDef"\'

echo -n Pushing MWClient defs, this rig to $1:
~/ExperimentXML-git/scripts/experimentSaveDefaults.py
echo -n saved, 
scp -pq ~/Desktop/org.histed.MWClientSavedVars.plist $1:/tmp
ssh $1 ExperimentXML-git/scripts/experimentLoadDefaultsFromFile.py /tmp/org.histed.MWClientSavedVars.plist
echo done.

