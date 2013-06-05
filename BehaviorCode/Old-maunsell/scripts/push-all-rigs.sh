#!/bin/sh
# push-all-rigs.sh [pushDefaults]
#   pushDefaults: {1|0} default true: whether to copy MWClient defaults from rig1
#
# first arg is MaunsellMouseN.local


scriptDir="$HOME/ExperimentXML-git/scripts"
for a in MaunsellMouse2.local MaunsellMouse3.local MaunsellMouse4.local; do

    if [ "$1" != "0" ]; then
	$scriptDir/push-defaults-one-rig.sh "$a"
    fi

    # MH 121209: with advent of multiple expts we shouldn't do this automatically, use experimentLoadDefaultsFromFile
    $scriptDir/push-files-one-rig.sh "$a"
done
