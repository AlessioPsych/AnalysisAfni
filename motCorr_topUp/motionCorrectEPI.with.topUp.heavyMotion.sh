#!/bin/bash

EPIDIR=$1
TOPOUTCOMEDIR=$2
EPINUM=$3

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data'
echo 'Inputs:'
echo 'EPIDIR=$1, input directory, where the EPI lives'
echo 'TOPOUTCOMEDIR=$2, input directory, where the TOPUP outcome files lives'
echo 'EPINUM=$3, which EPI run to cogerister everything to'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/motCorr_topUp/motionCorrectEPI.with.topUp.heavyMotion.R $EPIDIR $TOPOUTCOMEDIR $EPINUM


