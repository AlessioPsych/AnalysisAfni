#!/bin/bash

INPUT=$1
MINPROPORTION=$2
OUTPUTFILE=$3

if [ -z "$1" ]
then
echo 'Removes clusters of isolated voxels, smaller than the MINPROPORTION parameter, useful to clean up segmentations'
echo 'Inputs:'
echo 'INPUT, input anatomy'
echo 'MINPROPORTIO, path to where the AFNI atlases live'
echo 'OUTPUTFILE, output filename.'
echo 'example call: removeSmallClusters.sh input.nii.gz 0.05 cleanOutput.nii.gz'
exit 1
fi
 
Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/removeSmallClusters.R \
$INPUT $MINPROPORTION $OUTPUTFILE