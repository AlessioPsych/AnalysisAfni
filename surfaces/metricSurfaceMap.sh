#!/usr/bin/env bash

METRIC=$1
SURFACEFOLDER=$2

if [ -z "$1" ]
then
echo 'computes surface metric on a surface set (directory)'
echo 'Inputs:'
echo 'METRIC=$1, afni SurfaceMetrics flag (with dash), see SurfaceMetrics help for specifics'
echo 'SURFACEFOLDER=$2, surface folder over which to compute the selected metric'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRSURFACES/metricSurfaceMap.R $METRIC $SURFACEFOLDER
