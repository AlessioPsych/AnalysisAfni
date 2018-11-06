#!/bin/bash

INPUTDIR=$1
TARGETVOLUME=$2
COREGMAT=$3
OUTPUTDIR=$4


if [ -z "$1" ]
then
echo 'Bash script to apply a coregistration matrix to nifti files in one dir,'
echo 'add the *coreg* suffix, creates the dir OUTPUTDIR and copies the files'
echo 'make sense to use this for time series files for single volumes just use the'
echo 'plain 3dAllineate'
echo 'INPUTDIR=$1, input directory where the nifti lives' 
echo 'TARGETVOLUME=$2, target volume for the coregistration' 
echo 'COREGMAT=$3, affine coregistration matrix that projects the epi files into the target volume' 
echo 'OUTPUTDIR=$2, output directory where the coregistered niftis are stored'
echo
echo 'example call: applyCoregistrationToEpi.sh nifti_deob/ coregistration/meanEpiSingleShot.nii.gz coregistration/tMat.1D nifti_coreg/'
echo 'remember to put the backslash next to the dir names, as shown above'
echo
echo 'note that the code autobox the target volume to limit the file size of the coregistered time series'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/coregistration/applyCoregistrationToEpi.R \
$INPUTDIR \
$TARGETVOLUME \
$COREGMAT \
$OUTPUTDIR

