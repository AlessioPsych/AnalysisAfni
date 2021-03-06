#!/bin/bash

MASKFILE=$1
TSFILE=$2
OUTPUTNAME=$3
SURROUND=$4
POLORT=$5
SAMPLINGTIME=$6
STIMTYPE=$7
FITTEDPARAMETERS=$8
PREDICTEDTS=$9


if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'MASKFILE=$1, mask file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
echo 'FLAGSURROUND=$4, run surround fit? 1=yes 0=no'
echo 'POLORT=$5, -polort parameter for detrending'
echo 'SAMPLINGTIME=$6, set it to 0.166'
echo 'STIMTYPE=$7: 1:eye moving; 2:eye fix; 3:eye fix border; 4:prf; 5:eye fix border disappear'
echo 'FITTEDPARAMETERS=$8: already fitted parameters'
echo 'PREDICTEDTS=$9: obtained ts'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/twoDimensionalFit_Glasgow_eyeFix_gain_parallel_twoSteps.R $MASKFILE $TSFILE $OUTPUTNAME $SURROUND $POLORT $SAMPLINGTIME $STIMTYPE $FITTEDPARAMETERS $PREDICTEDTS
