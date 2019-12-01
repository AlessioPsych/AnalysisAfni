#!/bin/bash

INPUTDIR=$1

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory, where the EPI lives'
exit 1
fi

afni_proc.py -subj_id motionCorrect \
	-dsets $INPUTDIR*.nii	\
	-blocks despike volreg

tcsh -xef proc.motionCorrect |& tee output.proc.motionCorrect
