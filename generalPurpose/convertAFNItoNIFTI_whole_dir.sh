#!/bin/bash

INPUTDIR=$1
OUTPUTDIR=$2

if [ -z "$1" ]
then
echo 'converts the motion corrected volreg+orig files in a directory into nifti' 
echo 'renames the files into 01.nii 02.nii etc'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory, where the volreg+orig files lives'
echo 'OUTPUTDIR=$2, output directory, where the nii files will live' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/convertAFNItoNIFTI_whole_dir.R $INPUTDIR $OUTPUTDIR

