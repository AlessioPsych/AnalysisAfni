#!/bin/bash

ROI=$1
NODEFACTR=$2
KERNELRADIUS=$3
SURFACENAME=$4
ROTATION=$5
OUTNAME=$6

if [ -z "$1" ]
then
echo
echo
echo 'computes the orthogonal axis representing the flattened mesh'
echo
echo
echo 'Inputs:'
echo 'ROI=$1, e.g. 1D roi '
echo 'NODEFACTR=$2, downsampling factor for the number of nodes, try with 7 to begin with, small number -> longer computation time'
echo 'KERNELRADIUS=$3, radius of the interpolating kernel over the surface, in mm'
echo 'SURFACENAME=$4 surface name, e.g. std.141.lh.smoothwm.gii'
echo 'ROTATION=$5 parameter to rotate map'
echo 'OUTNAME=$6 output filename, e.g. 5_0.2_axis_lh'
echo
echo 'example call: compSpectralDecompRoi_manuscriptFigure.sh testRoi.1D.roi 1 1.8 std.60.rh.smoothwm.gii 0 testDecomposition.rh'
echo
echo
exit 1
fi


Rscript $AFNI_TOOLBOXDIR/devel/compSpectralDecompRoi_manuscriptFigure.R \
$ROI \
$NODEFACTR \
$KERNELRADIUS \
$AFNI_INSTALLDIR \
$AFNI_TOOLBOXDIRSURFACES \
$SURFACENAME \
$ROTATION \
$OUTNAME






