#!/bin/bash

INPUT=$1
REFPHASE=$2

if [ -z "$1" ]
then
echo 'Bash script to compute phase specified maps from travelling wave analysis.'
echo 'INPUT=$1, travelling wave input, see file compTravellingWave.sh' 
echo 'REFPHASE=$2, reference phase' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/fitEpiData/compPhaseSpecifiedMap.R \
$INPUT $REFPHASE $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES
