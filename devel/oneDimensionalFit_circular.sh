#!/bin/bash

EPIFILE=$1
TSFILE=$2
OUTPUTNAME=$3

if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'EPIFILE=$1, mean epi file'
echo 'TSFILE=$2, detrended ts file' 
echo 'OUTPUTNAME=$3, output name'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/oneDimensionalFit.R $EPIFILE $TSFILE $OUTPUTNAME $AFNI_INSTALLDIR 
