#!/bin/bash

INPUTEPIDIR=$1
INPUTEPIPHASEDIR=$2
EPIFORMAT=$3
NCORES=$4
OUTDIRNAME=$5

if [ -z "$1" ]
then
echo 'data denoising based on Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406, doi: 10.1016/j.neuroimage.2016.08.016'
echo 'it requires mrtrix3 installed, possibly via miniconda'
echo 'remember to - conda activate - before use'
echo 'Inputs:'
echo 'INPUTEPIDIR=$1, input EPI directory, put the backslash!'
echo 'INPUTEPIPHASEDIR=$2, input EPI PHASE directory, put the backslash!'
echo 'EPIFORMAT=$3, EPI format files, for example: *.nii'
echo 'NCORES=$4, number of cores to use' 
echo 'OUTDIRNAME=$5, output directory name'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/denoiseData_with_phase.R $INPUTEPIDIR $INPUTEPIPHASEDIR $EPIFORMAT $NCORES $OUTDIRNAME
