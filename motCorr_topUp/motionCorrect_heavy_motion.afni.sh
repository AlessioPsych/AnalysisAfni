#!/bin/bash

INPUTDIR=$1

if [ -z "$1" ]
then
echo 'computes despike and motion corrected epi data'
echo 'to use if data has a lot of motion'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory, where the EPI lives'
exit 1
fi


afni_proc.py -subj_id motionCorrect \
	-dsets EPI/*.nii \
	-blocks despike volreg \
	-volreg_opts_vr -twopass \
	-volreg_opts_vr -cubic \
	-volreg_opts_vr -maxite 38 \
	-volreg_opts_vr -final heptic \
	-volreg_interp -Fourier


