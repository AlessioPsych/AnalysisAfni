#!/usr/bin/env bash

DIR=$1
SMOOTHING=$2

if [ -z "$1" ]
 then
  echo 'Inputs:'
  echo 'DIR=$1'
  echo 'SMOOTHING=$2'
  exit 1
fi

Rscript $AFNI_TOOLBOXDIRSURFACES/smoothSurfaceMap.R $DIR $SMOOTHING
