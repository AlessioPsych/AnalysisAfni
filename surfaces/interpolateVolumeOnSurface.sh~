#!/bin/bash

if [ -z "$1" ]
 then
  echo 'link to bash script to interpolate a 3d volume on a surface, the volume'
  echo 'can 4D; interpolation is NN, you must have the file spec.surfaces.smoothed'
  echo 'in your folder as well as the boundaries files. Check the file defineBoundaries.sh'
  echo 'and generateSurfaces.sh for that. The outcome is a folder, each file in the folder'
  echo 'represents the interpolated values for a give boundary. Inputs:'
  echo 'VOLUME=$1, volume filename' 
  exit 1
fi

VOLUME=$1

Rscript $AFNI_TOOLBOXDIRSURFACES/interpolateSurface.R $VOLUME
