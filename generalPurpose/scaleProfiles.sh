#!/bin/bash

INPUTVOL=$1
NEWMAX=$2
NEWMIN=$3
OUTPUTFILENAME=$4

if [ -z "$1" ]
then
echo 'Bash script to rescale a profile volume generated from nighres or the MIPAV pipeline,'
echo 'the scaling is applied only to valid voxels (where he profile is not equal to zero). Inputs:'
echo 'INPUTVOL=$1, volume filename' 
echo 'NEWMAX=$2'
echo 'NEWMIN=$3'
echo 'OUTPUTFILENAME=$4, output volume filename' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/scaleProfiles.R \
 $INPUTVOL $NEWMAX $NEWMIN $OUTPUTFILENAME $AFNI_TOOLBOXDIRGENERALPURPOSE $AFNI_INSTALLDIR
