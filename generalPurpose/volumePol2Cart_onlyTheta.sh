#!/bin/bash

INPUTVOLUME=$1
THETAIDX=$2
OUTNAME=$3


if [ -z "$1" ]
then
echo 'converts polar coordinates (only theta) stored in a 4D volume (theta only) into cartesian coordinates. In this case rho is set to 1 by default to allow one dimensional conversion. N.B. in this case theta is given in degrees, not radiants. This is confusing, but is needed for the specific application for which this code has been developed.',
echo 'Inputs:'
echo 'INPUTVOLUME=$1, input volume, probably prf estimates or kinematic estimates'
echo 'THETAIDX=$2, volume index in afni where theta position (polar, is stored), afni starts counting from 0' 
echo 'OUTNAME=$3, outputname (with extension: e.g.: .nii.gz)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/volumePol2Cart_onlyTheta.R $INPUTVOLUME $THETAIDX $OUTNAME
