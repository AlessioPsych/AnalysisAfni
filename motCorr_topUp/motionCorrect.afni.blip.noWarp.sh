#!/bin/bash

EPIDIR=$1
TOPDIR=$2
EPINUM=$3
TOPNUM=$4
EPILIST=$5
TOPLIST=$6

if [ -z "$1" ]
then
echo 'computes the basic files for top up correction, blip up and blip down, to feed into afni.proc.py, it does NOT perform non linear warping, afni does that later'
echo 'Inputs:'
echo 'EPIDIR=$1, input directory, where the EPI lives'
echo 'TOPDIR=$2, input directory, where the TOPUP files lives'
echo 'EPINUM=$3, which EPI to coregister everything to?'
echo 'TOPNUM=$4, which EPI to coregister everything to?'
echo 'EPILIST=$5, which EPIS to compute the top up? e.g 1-3-7'
echo 'TOPLIST=$6, which TOP-UPS to compute the top up? e.g 1-3-7'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/motCorr_topUp/motionCorrect.afni.blip.noWarp.R $EPIDIR $TOPDIR $EPINUM $TOPNUM $EPILIST $TOPLIST


