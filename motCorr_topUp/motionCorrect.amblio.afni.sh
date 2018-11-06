#!/bin/bash

INPUTDIR=$1

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data for amblio dataset, removing 8 trs'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory, where the EPI lives'
exit 1
fi

afni_proc.py -subj_id motionCorrect \
	-dsets $INPUTDIR*.nii	\
	-blocks tcat despike volreg -tcat_remove_first_trs 8 -volreg_interp -Fourier 

tcsh -xef proc.motionCorrect |& tee output.proc.motionCorrect
