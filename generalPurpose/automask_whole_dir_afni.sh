#!/bin/bash

FILESDIR=$1
CLFRAC=$2

if [ -z "$1" ]
then
echo 'computes 3dAutomask on *.nii files in a directory, renames the files into 01.nii 02.nii etc etc'
echo 'Inputs:'
echo 'FILESDIR=$1, input directory, where the *.nii files lives'
echo 'CLFRAC=$2, fractional intensity threshold (0->1); smaller values give larger brain estimate' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/automask_whole_dir_afni.R $FILESDIR $CLFRAC

