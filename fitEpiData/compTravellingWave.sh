#!/bin/bash

INPUTTS=$1
CYCLES=$2
DETPOLORT=$3
PHASESHIFT=$4

if [ -z "$1" ]
then
echo 'Bash script to compute the travelling wave analysis.'
echo 'INPUTTS=$1, volume filename, time serie' 
echo 'CYCLES=$2, for how many cycles was the stimuli presented?'
echo 'DETPOLORT=$3, polort parameter in 3dDetrend'  
echo 'PHASESHIFT=$4, phase shift, in radiants, to match the beginning of the stimulus'
echo 'the code detrends the ts and compute the %signal change before fft'
echo 'as output you get a file trWave.nii.gz with fields: amplitude, coherence and phase'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/fitEpiData/compTravellingWave.R \
$INPUTTS $CYCLES $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES $DETPOLORT $PHASESHIFT $AFNI_TOOLBOXDIRGENERALPURPOSE

