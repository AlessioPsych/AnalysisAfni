#!/bin/bash

maindir="/home/fracasso/Data/openNeuro/ds001751_FlankerContextualModulation"
targetdir="derivatives/mrtrix3"

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

	cd "$maindir"

	echo "Processing directory: $dir"
	
#	if [ -d $maindir/$targetdir/$dir ]; then
#	        echo "removing folder $maindir/$targetdir/$dir" 
#		rm -R $maindir/$targetdir/$dir
#	fi
#	mkdir $maindir/$targetdir/$dir
	
	echo "data folder: $maindir/$dir/func"
	echo "target folder: $maindir/$targetdir/$dir/func"
	
	# copy func folder to denoise
#	echo "cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func"	
#	cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func
	
	cd "$maindir/$targetdir/$dir/func"
	echo "Current folder: $PWD"

	if [ -f ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz ]; then 
		echo "rm ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz"
		rm ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz 
	fi
	if [ -f ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz ]; then 
		echo "rm ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz"
		rm ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz 
	fi

	echo "3dWarp -NN -deoblique -prefix ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz"
	echo "3dWarp -NN -deoblique -prefix ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz"

	3dWarp -NN -deoblique -prefix ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect_deob.nii.gz ${dir}_task-flanker_run-1_bold_denoised_tSliceCorrect.nii.gz	
	3dWarp -NN -deoblique -prefix ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect_deob.nii.gz ${dir}_task-flanker_run-2_bold_denoised_tSliceCorrect.nii.gz

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


