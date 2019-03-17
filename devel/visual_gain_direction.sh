#!/bin/bash

MASKFILE=$1
TSFILE=$2
OUTPUTNAME=$3
FLAGSURROUND=$4
POLORT=$5
SAMPLINGTIME=$6
FITTYPE=$7

if [ -z "$1" ]
then
echo 'computes fit with visual + gain + direction parameters'
echo 'Inputs:'
echo 'MASKFILE=$1, mask file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
echo 'FLAGSURROUND=$4, run fine or coarse fit ? 1=fine 0=coarse'
echo 'POLORT=$5, -polort parameter for detrending'
echo 'SAMPLINGTIME=$6, set it to 0.166'
echo 'FITTYPE=$7, 1 (gain); 2 (only prf)'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim) library( matrixcalc )'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/visual_gain_direction.R $MASKFILE $TSFILE $OUTPUTNAME $FLAGSURROUND $POLORT $SAMPLINGTIME $FITTYPE
