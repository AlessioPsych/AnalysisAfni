#!/bin/bash

if [ -z "$1" ]
 then
  echo 'SURFACESDIR=$1, surfaces directory' 
  echo 'VOLUME=$2, volume filename' 
  exit 1
fi

SURFACESDIR=$1
VOLUME=$2


Rscript $AFNI_TOOLBOXDIRSURFACES/afniSurface.R $SURFACESDIR $VOLUME

##!/usr/bin/env bash

#MPRAGE=$1

#afni -niml &
#suma -spec surfaces_folder/spec.surfaces.smoothed -sv $MPRAGE &


