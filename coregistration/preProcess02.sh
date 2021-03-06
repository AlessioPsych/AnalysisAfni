#!/bin/bash

INPUTVOL=$1
PERCENTILE=$2
CLFRAC=$3
DILATE=$4
ERODE=$5
ZEROPAD=$6
OUTVOLUME=$7
UNIFTYPE=$8
KEEPOBLIQUE=$9

if [ -z "$1" ]
 then
  echo 'Bash script to pre process volume a volume before coregistration. Inputs:'
  echo 'INPUTVOL=$1, volume filename' 
  echo 'PERCENTILE=$2, masking based on intensity percentile [0 100]'
  echo 'CLFRAC=$3, 3dAutomask parameter [0, 1]'
  echo 'DILATE=$4, 3dAutomask parameter'
  echo 'ERODE=$5 3dAutomask paramete'
  echo 'ZEROPAD=$6, Zeropad, n slices, all directions'
  echo 'OUTVOLUME=$7, output volume'
  echo 'UNIFTYPE=$8, make volume uniform [0 (no), 1 (yes, 3dUnifize), 2 (yes, 3dUnifize inverted contrast), 3 (yes, 3dUniformize) ]'
  echo 'KEEPOBLIQUE=$9, keep oblique information in header [0 1], if you save nii or nii.gz, the header information is lost anyway'
  exit 1
fi


Rscript $AFNI_TOOLBOXDIR/coregistration/preProcess02.R \
 $INPUTVOL \
 $PERCENTILE \
 $CLFRAC \
 $DILATE \
 $ERODE \
 $ZEROPAD \
 $OUTVOLUME \
 $UNIFTYPE \
 $KEEPOBLIQUE

rm temp*.BRIK
rm temp*.HEAD
rm temp*.1D
