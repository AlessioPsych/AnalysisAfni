#!/bin/bash

INPUTANAT=$1
ATLASPATH=$2
NONLINFLAG=$3

if [ -z "$1" ]
then
echo 'Removes the cerebellum and brainstem from input anatomy, generates left hemisphere mask as well'
echo 'Inputs:'
echo 'INPUTANAT=$1, input anatomy'
echo 'ATLASPATH=$2, path to where the AFNI atlases live'
echo 'NONLINFLAG=$3, flag to enable non-linear coregistration (generally a good idea).'
echo 'example call: removeCerebellum.sh anatomy.nii.gz /packages/afni/17.0.13 1'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/removeCerebellum.R $INPUTANAT $ATLASPATH $NONLINFLAG

