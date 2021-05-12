#!/bin/bash
# bash file to call segment_white_matter.R via Rscript

ANATFILE=$1
WHITEMATTERMASK=$2
THR=$3

if [ -z "$1" ]
then
echo 'Bash script to segment T1 white matter. Inputs:'
echo 'ANATFILE=$1, matrix filename' 
echo 'WHITEMATTERMASK=$2, white matter mask generated from segment_white_matter.sh' 
echo 'THR=$3, threshold for gm boundary, four values in this order:' 
echo '1=frontal, 2=temporal, 3=parietal, 4=occipital,'
echo 'larger values give thicker estimates, for example 0.65-0.65-0.55-0.65'
echo
exit 1
fi

Rscript $AFNI_TOOLBOXDIRCOREGISTRATION/segment_gray_matter_two_steps.R $ANATFILE $WHITEMATTERMASK $AFNI_INSTALLDIR $AFNI_TOOLBOXDIRSURFACES $AFNI_ATLASDIR_TOOLBOX $THR
