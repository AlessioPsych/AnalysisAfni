#!/bin/bash

ROIFILE=$1
PROFILEFILE=$2
OUTPUTNAME=$3


if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'ROIFILE=$1, roi file file, selected voxels, probably something thresholded for cortical depth, to avoid multiple estimation of similar profiles'
echo 'PROFILEFILE=$2, profiles to model, probably from volumetric layering' 
echo 'OUTPUTNAME=$3, output name, without the extension, something like delMe'
echo 'it needs the following R libraries:'
echo 'library( parallel )'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/modelProfiles.R $ROIFILE $PROFILEFILE $OUTPUTNAME
