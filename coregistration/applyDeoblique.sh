#!/bin/bash

INPUTDIR=$1
OUTPUTDIR=$2

if [ -z "$1" ]
then
echo 'Bash script to deoblique all the nifti files in one dir,'
echo 'add the *deob* suffix, creates the dir OUTPUTDIR and copies the files'
echo 'INPUTDIR=$1, input directory where the nifti lives' 
echo 'OUTPUTDIR=$2, output directory where the deoblique niftis are stored'
echo
echo 'example call: applyDeoblique.sh niftiFiles_delme/ nifti_deob/'
echo 'remember to put the backslash next to the dir names, as shown above'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/coregistration/applyDeoblique.R \
$INPUTDIR \
$OUTPUTDIR