#!/bin/bash

INPUTVOLUME=$1
DEPTH=$2
OUTNAME=$3

if [ -z "$1" ]
then
echo
echo
echo 'this is very specific for an application to anatomical data and hemianopia, for visualization'
echo 'purposes only, it expands volumetric data of a surface, along cortical depth, do not use for'
echo 'computation, just visualization'
echo
echo 'Inputs:'
echo 'INPUTVOLUME=$1, values on a surface in volumetric space, to expand across depth'
echo 'DEPTH=$2, e.g. depth file, in the same master as the input, as usuals'
echo 'OUTNAME=$7, file output name'
echo
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/surfaces/expandVolumeValuesAcrossDepth.R \
$INPUTVOLUME \
$DEPTH \
$OUTNAME \
$AFNI_INSTALLDIR \
$AFNI_TOOLBOXDIRSURFACES
