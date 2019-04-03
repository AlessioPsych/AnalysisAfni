#!/bin/bash

ROIFILE=$1
OUTPUTNAME=$2
SCALING=$3

if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'ROIFILE=$1, roi file file, selected voxels,'
echo 'OUTPUTNAME=$2, output name, without the extension, something like delMe'
echo 'SCALING=$3, scaling factor for eccentricity map'
echo 'it needs the following R libraries:'
echo 'library( parallel )'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/convertAxisToMap.R $ROIFILE $OUTPUTNAME $SCALING
