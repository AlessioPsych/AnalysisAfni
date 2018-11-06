#!/bin/bash

INPUTVOL=$1
DETRENDPOLORT=$2
STIMULI1D=$3

if [ -z "$1" ]
then
echo 'Bash script to compute the phisiological noise predictors based on compCor.'
echo 'INPUTVOL=$1, volume filename' 
echo 'DETRENDPOLORT=$2, parameter to detrend the data in 3dDetrend, same as polort for AFNI' 
echo 'STIMULI1D=$3, .1D file with 1 column and as many rows as the number of trs in the INPUTVOL'
echo 'with 1 when the stimuli was presented and 0 when it was not'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/compCor.R \
$INPUTVOL $DETRENDPOLORT $STIMULI1D $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES
