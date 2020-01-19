#!/bin/bash

TRWAVEINPUT=$1
CYCLESHIFT=$2
TRWAVEOUTPUT=$3

if [ -z "$1" ]
then
echo 'Bash script to shift the phase from travelling wave analysis to account for hemodynamicl delay.'
echo 'TRWAVEINPUT=$1, travelling wave output volume filename, from "compTravellingWave.sh"' 
echo 'CYCLESHIFT=$2, shift, in degrees, not radiants'
echo 'TRWAVEOUTPUT=$3, output name'
echo 'example call: compShiftPhase.sh trWave.nii.gz 30 trWaveShift30.nii.gz'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/fitEpiData/compShiftPhase.R \
$TRWAVEINPUT $CYCLESHIFT $TRWAVEOUTPUT
