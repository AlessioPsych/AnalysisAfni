#!/bin/bash

INPUTEPIDIR=$1
EPIFORMAT=$2
INPUTJSONDIR=$3
OUTDIRNAME=$4

if [ -z "$1" ]
then
echo 'performs time slice correction based on the provided json files'
echo 'Inputs:'
echo 'INPUTEPIDIR=$1, input EPI directory, put the backslash!'
echo 'EPIFORMAT=$2, EPI format files, for example: *.nii'
echo 'INPUTJSONDIR=$3, input JSON directory, put the backslash!' 
echo 'OUTDIRNAME=$4, output directory name'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/timeSliceCorrection.R $INPUTEPIDIR $EPIFORMAT $INPUTJSONDIR $OUTDIRNAME
