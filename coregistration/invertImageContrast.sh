#!/bin/bash

INPUTVOL=$1
OUTPUTFILENAME=$2

if [ -z "$1" ]
then
echo 'Bash script to invert the contrast of a give image. Inputs:'
echo 'INPUTVOL=$1, volume filename' 
echo 'OUTPUTFILENAME=$2, output volume filename' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/coregistration/invertImageContrast.R \
 $INPUTVOL \
 $OUTPUTFILENAME