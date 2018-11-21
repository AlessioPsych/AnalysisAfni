#!/bin/bash

SURFACEROIFILENAME=$1
ROIOUTPUTNAME=$2


if [ -z "$1" ]
then
echo 'converts an roi generated from 3dVol2Surf (.dset file with 7 columns)'
echo 'into a .1D file with a list of nodes (column 1) and a list of ones (column 2)'
echo 'useful to pass from volume based rois ro surface based rois for, for example'
echo 'axis computation with compSpectralDecomposition.sh'
echo 'it takes the last column of .dset file, finds the ones and reports the corresponding'
echo 'first column (roi nodes)'
echo 'Inputs:'
echo 'SURFACEROIFILENAME=$1, surface filename (typically .dset, 7 columns)'
echo 'ROIOUTPUTNAME=$2, roi output filename' 
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/sur2vol_2_1dROI.R $SURFACEROIFILENAME $ROIOUTPUTNAME
