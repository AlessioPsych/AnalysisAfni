#!/bin/bash

ROI=$1
EPI=$2
CMMATRIX=$3
ALMATRIX=$4
NAME=$5
INTERPMETHOD=$6

roiFilename=$(printf '%s_al_epi+orig' ${NAME} )

3dcopy ${ROI} ROI+orig
3dcopy ${EPI} EPI+orig

3dZeropad -z 10 -prefix ROI_zp+orig ROI+orig
3dZeropad -z 10 -prefix EPI_zp+orig EPI+orig

3dAllineate   -1Dmatrix_apply $CMMATRIX    \
              -prefix MPRAGE_zp_cm+orig -master EPI_zp+orig \
              -input ROI_zp+orig \
	      -final $INTERPMETHOD

3dAllineate -final $INTERPMETHOD -1Dmatrix_apply $ALMATRIX -prefix ROI_no_zp+orig MPRAGE_zp_cm+orig

3dZeropad -z -10 -prefix ${roiFilename} ROI_no_zp+orig

3dAFNItoNIFTI ${roiFilename}
gzip $(printf '%s_al_epi.nii' ${NAME} )

rm MPRAGE*.BRIK
rm MPRAGE*.HEAD
rm *al_epi*.BRIK
rm *al_epi*.HEAD
rm EPI*.BRIK
rm EPI*.HEAD
rm ROI*.BRIK
rm ROI*.HEAD


