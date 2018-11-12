#!/bin/bash

if [ -z "$1" ]
 then
  echo 'start from an initialized set of surfaces (boundaries files,'
  echo 'check the file defineBoundaries.sh and generateSurfaces.sh for that), grows a normal'
  echo 'from the selected surface to the CSF boundary and to the GM/WM boundary, for each node in the surface, and'
  echo 'stores the mapping in a separate folder: Inputs:'
  echo 'STARTINGSURFACE=$1' 
  exit 1
fi

STARTINGSURFACE=$1

Rscript $AFNI_TOOLBOXDIRSURFACES/mappingBetweenSurfaces_WHOLE.R $STARTINGSURFACE
