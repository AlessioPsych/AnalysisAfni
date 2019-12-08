#!/bin/bash

ROIFILE=$1
PROFILEFILE=$2
OUTPUTNAME=$3
POLYPAR=$4


if [ -z "$1" ]
then
echo 'computes 1d fit'
echo 'Inputs:'
echo 'ROIFILE=$1, roi file file, selected voxels, probably something thresholded for cortical depth, to avoid multiple estimation of similar profiles'
echo 'PROFILEFILE=$2, profiles to model, probably from volumetric layering' 
echo 'OUTPUTNAME=$3, output name, without the extension, something like delMe'
echo 'POLYPAR=$4, fitting polinomial parameter along cortical depth (1: linear, 2: quadratic etc etc)'
echo 'it needs the following R libraries:'
echo 'library( parallel )'
echo 'NB: profiles are scaled within the ROI between 1000 and 8000, scaled profiles are saved'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/modelProfiles_onlyPoly.R $ROIFILE $PROFILEFILE $OUTPUTNAME $POLYPAR
