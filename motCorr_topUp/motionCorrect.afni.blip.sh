#!/bin/bash

EPIDIR=$1
TOPDIR=$2
EPINUM=$3
TOPNUM=$4
EPILIST=$5
TOPLIST=$6
MAXLEV=$7
BLUR=$8

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data'
echo 'Inputs:'
echo 'EPIDIR=$1, input directory, where the EPI lives'
echo 'TOPDIR=$2, input directory, where the TOPUP files lives'
echo 'EPINUM=$3, which EPI to coregister everything to?'
echo 'TOPNUM=$4, which EPI to coregister everything to?'
echo 'EPILIST=$5, which EPIS to compute the top up? e.g 1-3-7'
echo 'TOPLIST=$6, which TOP-UPS to compute the top up? e.g 1-3-7'
echo 'MAXLEV=$7, which minimum grid size for 3dQWarp (mm)? try with 5 to begin with'
echo 'BLUR=$8, how much to blur the data for 3dQwarp? use -1 in general (no blur), try to blur only for problematic volumes' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/motCorr_topUp/motionCorrect.afni.blip.R $EPIDIR $TOPDIR $EPINUM $TOPNUM $EPILIST $TOPLIST $MAXLEV $BLUR


