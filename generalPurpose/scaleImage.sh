#!/bin/bash

INPUTVOL=$1
NEWMAX=$2
NEWMIN=$3
OUTPUTFILENAME=$4
ROI=$5

if [ -z "$1" ]
then
echo 'Bash script to rescale a given image. Inputs:'
echo 'INPUTVOL=$1, volume filename' 
echo 'NEWMAX=$2'
echo 'NEWMIN=$3'
echo 'OUTPUTFILENAME=$4, output volume filename' 
echo 'ROI=$5, roi that delimits the scaling'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/scaleImage.R \
 $INPUTVOL $NEWMAX $NEWMIN $OUTPUTFILENAME $AFNI_TOOLBOXDIRGENERALPURPOSE $AFNI_INSTALLDIR $ROI
