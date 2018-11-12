#!/bin/bash

T13DEPI=$1
T13DEPI_TOPUP=$2
T13DEPI_MASK=$3
T13DEPI_TOPUP_MASK=$4
MINPATCH=$5

if [ -z "$1" ]
then
echo 'Inputs:'
echo 'T13DEPI=$1, imago to correct'
echo 'T13DEPI_TOPUP=$2, top up'
echo 'T13DEPI_MASK=$3, binary mask for T13DEPI, use 3dUnifize and 3dmask_tool to generate it'
echo 'T13DEPI_TOPUP_MASK=$4, binary mask for T13DEPI, use 3dUnifize and 3dmask_tool to generate it'
echo 'MINPATCH=$5, minimum patch size to optimize, in mm. E.g. 5'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/motCorr_topUp/motionCorrect.afni.T1wEPI.R $T13DEPI $T13DEPI_TOPUP $T13DEPI_MASK $T13DEPI_TOPUP_MASK $MINPATCH

