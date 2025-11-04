#!/bin/bash

# Instructions to run:
# sudo docker run -v /media/alessiofracasso/DATADRIVE1/Flanker:/mrtrixDataFolder --rm -it mrtrix3/mrtrix3
# cd /mrtrixDataFolder
# sh 01_00_denoiseData.sh

maindir="/mrtrixDataFolder"
targetdir="derivatives/mrtrix3"

echo "$maindir"
echo "$targetdir"

cd "$maindir"
echo "Current folder: $PWD"

# print list of participants
\ls -d sub-* > subdirs.txt

if [ -d $maindir/$targetdir ]; then
	echo "removing folder $maindir/$targetdir" 
	rm -R $maindir/$targetdir
fi
mkdir $maindir/$targetdir

while IFS= read -r dir; do

	cd "$maindir"

	echo "Processing directory: $dir"
	
	if [ -d $maindir/$targetdir/$dir ]; then
	        echo "removing folder $maindir/$targetdir/$dir" 
		rm -R $maindir/$targetdir/$dir
	fi
	mkdir $maindir/$targetdir/$dir
	
	echo "data folder: $maindir/$dir/func"
	echo "target folder: $maindir/$targetdir/$dir/func"
	
	# copy func folder to denoise
	echo "cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func"	
	cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func
	
	# denoise data
	cd "$maindir/$targetdir/$dir/func"
	echo "Current folder: $PWD"

	echo "dwidenoise -nthreads 4 ${dir}_task-flanker_run-1_bold.nii.gz ${dir}_task-flanker_run-1_bold_denoised.nii.gz"
	echo "dwidenoise -nthreads 4 ${dir}_task-flanker_run-2_bold.nii.gz ${dir}_task-flanker_run-2_bold_denoised.nii.gz"
	dwidenoise -nthreads 4 ${dir}_task-flanker_run-1_bold.nii.gz ${dir}_task-flanker_run-1_bold_denoised.nii.gz
	dwidenoise -nthreads 4 ${dir}_task-flanker_run-2_bold.nii.gz ${dir}_task-flanker_run-2_bold_denoised.nii.gz
	rm ${dir}_task-flanker_run-1_bold.nii.gz
	rm ${dir}_task-flanker_run-2_bold.nii.gz

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


