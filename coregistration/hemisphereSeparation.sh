#!/bin/bash

INPUTANAT=$1
INPUTSEG=$2
ATLASPATH=$3
NONLINFLAG=$4

if [ -z "$1" ]
then
echo 'Separates a segmentation into left and right hemisphere'
echo 'Inputs:'
echo 'INPUTANAT=$1, input anatomy'
echo 'INPUTSEG=$2, input segmentation'
echo 'ATLASPATH=$3, atlas to coregister the anatomy to'
echo 'NONLINFLAG=$4, flag to enable non-linear hemishpere separation.'
echo 'example call: hemisphereSeparation.sh anatomy.nii.gz segmentation.nii.gz /usr/share/afni/atlases/TT_icbm452+tlrc 1'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/hemisphereSeparation.R $INPUTANAT $INPUTSEG $ATLASPATH $NONLINFLAG

