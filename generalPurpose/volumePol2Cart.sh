#!/bin/bash

INPUTVOLUME=$1
THETAIDX=$2
RHOIDX=$3
OUTNAME=$4


if [ -z "$1" ]
then
echo 'converts polar coordinates stored in a 4D volume (theta and rho) into cartesian coordinates',
echo 'Inputs:'
echo 'INPUTVOLUME=$1, input volume, probably prf estimates or kinematic estimates'
echo 'THETAIDX=$2, volume index in afni where theta position (polar, is stored), afni starts counting from 0' 
echo 'RHOIDX=$3, volume index in afni where rho position (polar, is stored), afni starts counting from 0'
echo 'OUTNAME=$4, outputname (with extension: e.g.: .nii.gz)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/volumePol2Cart.R $INPUTVOLUME $THETAIDX $RHOIDX $OUTNAME
