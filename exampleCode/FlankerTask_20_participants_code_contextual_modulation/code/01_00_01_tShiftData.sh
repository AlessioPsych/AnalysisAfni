#!/bin/bash

maindir="/mnt/disk01/ds001751_FlankerTask_context/derivatives/mrtrix3"
#targetdir="derivatives/mrtrix3"

echo "$maindir"
echo "$targetdir"

cd "$maindir"
echo "Current folder: $PWD"

# print list of participants
\ls -d sub-* > subdirs.txt

#if [ -d $maindir/$targetdir ]; then
#	echo "removing folder $maindir/$targetdir" 
#	rm -R $maindir/$targetdir
#fi
#mkdir $maindir/$targetdir

while IFS= read -r dir; do

	cd "$maindir/$dir/func"

	echo "Processing directory: $dir"
	echo
	echo "folder: $PWD"
	echo
	echo "slice time correction for: ${dir}_task-flanker_run-1_bold_denoised.nii.gz"
	echo "slice time correction with: ${dir}_task-flanker_run-1_bold.json"
	echo
	echo "slice time correction for: ${dir}_task-flanker_run-2_bold_denoised.nii.gz"
	echo "slice time correction with: ${dir}_task-flanker_run-2_bold.json"

	echo "Rscript $AFNI_TOOLBOXDIRGENERALPURPOSE/timeSliceCorrection_singleEpi_noRJSON.R ${dir}_task-flanker_run-1_bold_denoised.nii.gz ${dir}_task-flanker_run-1_bold.json ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz"
	echo
	echo "Rscript $AFNI_TOOLBOXDIRGENERALPURPOSE/timeSliceCorrection_singleEpi_noRJSON.R ${dir}_task-flanker_run-2_bold_denoised.nii.gz ${dir}_task-flanker_run-2_bold.json ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz"

	Rscript $AFNI_TOOLBOXDIRGENERALPURPOSE/timeSliceCorrection_singleEpi_noRJSON.R ${dir}_task-flanker_run-1_bold_denoised.nii.gz ${dir}_task-flanker_run-1_bold.json ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz
	
	Rscript $AFNI_TOOLBOXDIRGENERALPURPOSE/timeSliceCorrection_singleEpi_noRJSON.R ${dir}_task-flanker_run-2_bold_denoised.nii.gz ${dir}_task-flanker_run-2_bold.json ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


