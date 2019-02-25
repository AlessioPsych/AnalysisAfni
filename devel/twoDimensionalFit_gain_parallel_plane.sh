#!/bin/bash

MASKFILE=$1
TSFILE=$2
OUTPUTNAME=$3
FLAGFINEFIT=$4
POLORT=$5
SAMPLINGTIME=$6
STIMTYPE=$7
FITSURROUND=$8

if [ -z "$1" ]
then
echo 'computes 2d prf fit on different stimuli'
echo 'Inputs:'
echo 'MASKFILE=$1, mask file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
echo 'FLAGFINEFIT=$4, run fine or coarse fit ? 1=fine 0=coarse'
echo 'POLORT=$5, -polort parameter for detrending'
echo 'SAMPLINGTIME=$6, set it to 0.166'
echo 'STIMTYPE=$7: 1:eye moving; 2:eye fix; 3:eye fix border; 4:prf; 5:eye fix border disappear'
echo 'FITSURROUND=$8, fit surround? 1=yes, 2=no'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/twoDimensionalFit_Glasgow_eyeFix_gain_parallel_plane.R $MASKFILE $TSFILE $OUTPUTNAME $FLAGFINEFIT $POLORT $SAMPLINGTIME $STIMTYPE $FITSURROUND
