#!/bin/bash

EPIFILE=$1
TSFILE=$2
OUTPUTNAME=$3
FLAGSURROUND=$4
POLORT=$5
FITINTERCEPT=$6

if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'EPIFILE=$1, mask file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
echo 'FLAGSURROUND=$4, run surround fit? 1=yes 0=no'
echo 'POLORT=$5, -polort parameter for detrending'
echo 'FITINTERCEPT=$6, fit intercept to detrended ts? probably you want this to 1'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/twoDimensionalFit.R $EPIFILE $TSFILE $OUTPUTNAME $AFNI_INSTALLDIR $FLAGSURROUND $POLORT $FITINTERCEPT
