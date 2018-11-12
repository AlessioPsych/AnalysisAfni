#!/bin/bash

INPUTVOL=$1
MASTERVOL=$2

if [ -z "$1" ]
then
echo 'Bash script to copy the oblique matrix from one dataset to another.'
echo 'be careful, this creates a file with a modified header, you have to know'
echo 'what you are doing, it copies the whole IJK_TO_DICOM_REAL matrix, so'
echo 'it is going to change also the center of the volume.'
echo 'INPUTVOL=$1, volume filename' 
echo 'MASTERVOL=$2, volume from which the matrix is going to be copied' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/coregistration/applyObliqueMatrix.R \
 $INPUTVOL \
 $MASTERVOL $AFNI_INSTALLDIR