#!/bin/bash

INPUTEPIDIR=$1
EPIFORMAT=$2
NCORES=$3
OUTDIRNAME=$4

if [ -z "$1" ]
then
echo 'data denoising based on Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406, doi: 10.1016/j.neuroimage.2016.08.016'
echo 'it requires mrtrix3 installed, possibly via miniconda'
echo 'remember to - conda activate - before use'
echo 'Inputs:'
echo 'INPUTEPIDIR=$1, input EPI directory, put the backslash!'
echo 'EPIFORMAT=$2, EPI format files, for example: *.nii'
echo 'NCORES=$3, number of cores to use' 
echo 'OUTDIRNAME=$4, output directory name'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/denoiseData.R $INPUTEPIDIR $EPIFORMAT $NCORES $OUTDIRNAME
