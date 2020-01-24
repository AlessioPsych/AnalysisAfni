#!/bin/bash

INPUT=$1
CYCLESHIFT=$2
OUTPUT=$3
VOLTOSHIFT=$4

if [ -z "$1" ]
then
echo 'Bash script to shift the phase from model analysis.'
echo 'INPUT=$1, travelling wave output volume filename, from "compTravellingWave.sh"' 
echo 'CYCLESHIFT=$2, shift, in degrees, not radiants'
echo 'OUTPUT=$3, output name'
echo 'VOLTOSHIFT=$4, which volume in the input file to shift? in R convention, so starting from 1, not afni convention'
echo 'example call: compShiftPhase_model.sh model.nii.gz 30 model30.nii.gz 6'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/fitEpiData/compShiftPhase_model.R \
$INPUT $CYCLESHIFT $OUTPUT $VOLTOSHIFT
