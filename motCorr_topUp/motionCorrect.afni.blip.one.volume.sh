#!/bin/bash

EPIDIR=$1
TOPDIR=$2
EPINUM=$3
TOPNUM=$4

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data'
echo 'Inputs:'
echo 'EPIDIR=$1, input directory, where the EPI lives'
echo 'TOPDIR=$2, input directory, where the TOPUP files lives'
echo 'EPINUM=$3, which EPI to coregister everything to?'
echo 'TOPNUM=$4, which TOP-UP to use'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/motCorr_topUp/motionCorrect.afni.blip_singleVolume.R $EPIDIR $TOPDIR $EPINUM $TOPNUM


