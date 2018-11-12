#!/bin/bash
# bash file to call segment_white_matter.R via Rscript

ANATFILE=$1
CLUSTPROPLOCAL=$2
CLUSTPROPGLOBAL=$3
NITER=$4
DISTANCETHR=$5
BIAS=$6
PARAMWM=$7
PARAMCLEAN=$8

if [ -z "$1" ]
then
echo 'Bash script to segment T1 white matter. Inputs:'
echo 'ANATFILE=$1, matrix filename' 
echo 'CLUSTPROPLOCAL=$2, smallest proportion of voxels allowed to form an isolated cluster, locally, use values around [0.01 0.1], check the outcome!'
echo 'CLUSTPROPGLOBAL=$3, smallest proportion of voxels allowed to form an isolated cluster, globally, use values around [0.01 0.1], check the outcome!'
echo 'smaller numbers mean more inclusive threshold -> more white matter, decrease the number if there are missing portion of white matter'
echo 'NITER=$4, number of iterations for afni 3dSeg, this is important, start from ~20, it can be anything between 10 and 200 and beyond, check the results visually'
echo 'DISTANCETHR=$5, distance in voxel from WM over which a voxels cannot be assigned to WM, must be arranged like this e.g 8-8-8-7-7-7'
echo 'BIAS=$6, distance in voxel to grow (positive number) or erode (negative number) the final wm mask, just a global adjustment, normally set it to 0'
exit 1
fi

if [ -z "$7" ]
then
PARAMWM=0
PARAMCLEAN=1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/segment_white_matter_15082017_02.R $ANATFILE $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES $AFNI_ATLASDIR $CLUSTPROPLOCAL $CLUSTPROPGLOBAL $NITER $DISTANCETHR $AFNI_TOOLBOXDIRCOREGISTRATION $BIAS $PARAMWM $PARAMCLEAN
