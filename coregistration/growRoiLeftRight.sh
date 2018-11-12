#!/bin/bash
# smooths the left_right separation from the atlas

VOLUME=$1
SMOOTH=$2

if [ -z "$1" ]
 then
  echo 'smooths and expands the left_right roi ported from the atlas. Inputs:'
  echo 'VOLUME=$1, left_right roi volume filename' 
  echo 'SMOOTH=$2, smoothing parameter used by 3dmerge (afni)'
  exit 1
fi

Rscript $AFNI_TOOLBOXDIR/coregistration/growRoiLeftRight.R $VOLUME $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES $SMOOTH
