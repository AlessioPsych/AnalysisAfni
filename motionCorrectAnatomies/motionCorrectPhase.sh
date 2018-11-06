#!/usr/bin/env bash

INDIR=$1
SEARCHSTRING=$2
MATSTRING=$3

if [ -z "$1" ]
then
echo 'motion corrects volumes based on provided matrices (from AFNI), recursively'
echo 'Inputs:'
echo 'INDIR=$1, directory where the volumes are stored'
echo 'SEARCHSTRING=$2, regular expr to search for specific volumes in the directory'
echo 'MATSTRING=$3, regular expr to search for specific matrix files in the directory'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/motionCorrectPhase.R $INDIR $SEARCHSTRING $MATSTRING
