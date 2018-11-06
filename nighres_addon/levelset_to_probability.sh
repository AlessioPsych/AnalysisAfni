#!/bin/bash

INPUTVOL=$1
OUTPUTFILENAME=$2
THR=$3

if [ -z "$1" ]
then
echo 'Bash script convert a nighres levelset to probability.' 
echo 'it requires the library pracma in R. Inputs:'
echo 'INPUTVOL=$1, volume filename' 
echo 'OUTPUTFILENAME=$2, output volume filename' 
echo 'THR=$3, sigmoid stepness'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/nighres_addon/levelset_to_probability.R \
$INPUTVOL $OUTPUTFILENAME $THR $AFNI_INSTALLDIR
