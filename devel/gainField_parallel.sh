#!/bin/bash

MASKFILE=$1
TSFILE=$2
OUTPUTNAME=$3
FLAGSURROUND=$4
POLORT=$5
FITINTERCEPT=$6
SAMPLINGTIME=$7
STIMTYPE=$8


if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'MASKFILE=$1, mask file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
echo 'FLAGSURROUND=$4, run surround fit? 1=yes 0=no'
echo 'POLORT=$5, -polort parameter for detrending'
echo 'FITINTERCEPT=$6, fit intercept to detrended ts? probably you want this to 1'
echo 'SAMPLINGTIME=$7, set it to 0.166'
echo 'STIMTYPE=$8: 1:eye moving; 2:eye fix; 3:eye fix border; 4:prf; 5:eye fix border disappear'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/gainField_parallel.R $MASKFILE $TSFILE $OUTPUTNAME $AFNI_INSTALLDIR $FLAGSURROUND $POLORT $FITINTERCEPT $SAMPLINGTIME $STIMTYPE 
