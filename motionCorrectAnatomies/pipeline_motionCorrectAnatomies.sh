#!/bin/bash

if [ -z "$1" ]
then
echo 'motion correct a set of MPRAGE anat and pd correct'
echo 'remember to create 2 directories:'
echo 'originalNifti and pdVolumes before running it'
echo 'example call: pipeline_motionCorrectAnatomies.sh 1'
exit 1
fi

motionCorrect.sh

computeMeanMotionCorrected.sh

preprocessingPd.sh
