#!/usr/bin/env bash

INPUTDIR=$1
SMOOTHING=$2
SURFACESDIR=$3

if [ -z "$1" ]
 then
  echo 'Inputs:'
  echo 'INPUTDIR=$1, dir where surfaces are stored'
  echo 'SMOOTHING=$2 smoothing kernel'
  echo 'SURFACESDIR=$3 , remember to put the backslash'
  exit 1
fi

Rscript $AFNI_TOOLBOXDIRSURFACES/smoothSurfaceMap.R $INPUTDIR $SMOOTHING $SURFACESDIR
