#!/bin/bash

INPUTVOLUME=$1
XIDX=$2
YIDX=$3
OUTNAME=$4


if [ -z "$1" ]
then
echo 'converts cartesian coordinates stored in a 4D volume (x and y) into polar coordinates',
echo 'useful for interpolating prf parameters in antomical grid and avoiding the problem of averaging'
echo 'neighboring voxels in polar coordinates'
echo 'Inputs:'
echo 'INPUTVOLUME=$1, input volume, probably prf estimates or kinematic estimates'
echo 'XIDX=$2, volume index in afni where x position (cartesian, is stored), afni starts counting from 0' 
echo 'YIDX=$3, volume index in afni where y position (cartesian, is stored), afni starts counting from 0'
echo 'OUTNAME=$4, outputname (with extension: e.g.: .nii.gz)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/volumeCart2Pol.R $INPUTVOLUME $XIDX $YIDX $OUTNAME
