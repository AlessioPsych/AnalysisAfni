#!/usr/bin/env bash

NIIFILEIN=$1
MPRAGE=$2
TMAT=$3

if [ -z "$1" ]
then
echo 'computes angle between B0 and normals to surface'
echo 'Inputs:'
echo 'NIIFILEIN=$1, input nifti file'
echo 'MPRAGE=$2, mprage to align centers to, possibly deoblique before running this'
echo 'TMAT=$3, t mat from freeview'
exit 1
fi


Rscript $AFNI_TOOLBOXDIR/generalPurpose/freeview_transformNiftyfromMat.R $NIIFILEIN $MPRAGE $TMAT



