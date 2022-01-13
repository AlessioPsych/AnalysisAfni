#!/bin/bash

ROTATION=$1
OUTNAME=$2
INPUTFILE=$3

if [ -z "$1" ]
then
echo
echo
echo 'rotates a spectral decomposition file'
echo
echo
echo 'Inputs:'
echo 'ROTATION=$1 parameter to rotate map, in radiants'
echo 'OUTNAME=$2 output filename, e.g. 5_0.2_axis_lh'
echo 'INPUTFILE=$3 input file, spectral decomposition results (nodes, first component, second component) '
echo
echo 'example call: rotSpectralDecompRoi_Coimbra.sh 0.5 testRoi 7_0.4_axis_lh_axis.dset'
echo
echo
exit 1
fi


Rscript $AFNI_TOOLBOXDIR/devel/rotSpectralDecompRoi_Coimbra.R \
$ROTATION \
$OUTNAME \
$INPUTFILE





