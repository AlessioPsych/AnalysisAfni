#!/bin/bash

INPUTMASK=$1
DISTANCETHR=$2
OUTNAME=$3

if [ -z "$1" ]
then
echo 'Erodes the mask, computes distance between resulting mask and eroded voxels, eliminates all the voxels exceeding the euclidian distance in DISTANCETHR'
echo 'useful to clean up the masks for segmentation purposes'
echo 'Inputs:'
echo 'INPUTMASK=$1, input mask'
echo 'DISTANCETHR=$2, threshold distance'
echo 'OUTNAME=$3, output filename'
echo 'example call: erodeMask_distance.sh brainMask.nii.gz 3 brainMask_thr.nii.gz'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/erodeMask_distance.R $INPUTMASK $DISTANCETHR $OUTNAME

