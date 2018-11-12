#!/usr/bin/env bash

METRIC=$1

if [ -z "$1" ]
then
echo 'computes surface metric on a surface set (directory)'
echo 'Inputs:'
echo 'METRIC=$1, afni SurfaceMetrics flag (with dash)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIRSURFACES/metricSurfaceMap.R $METRIC
