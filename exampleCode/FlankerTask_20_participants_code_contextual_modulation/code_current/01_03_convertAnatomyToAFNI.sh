#!/bin/bash

maindir="/home/fracasso/Data/openNeuro/ds001751_FlankerContextualModulation"
targetdir="derivatives/Freesurfer_output"

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
	
	#if [ -d $maindir/$targetdir/$dir ]; then
	#        echo "removing folder $maindir/$targetdir/$dir" 
	#		rm -R $maindir/$targetdir/$dir
	#fi
	#mkdir $maindir/$targetdir/$dir
	
	#echo "data folder: $maindir/$dir/func"
	#echo "target folder: $maindir/$targetdir/$dir/func"
	
	# copy func folder to denoise
	#echo "cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func"	
	#cp -R $maindir/$dir/func $maindir/$targetdir/$dir/func
	
	# denoise data
	cd "$maindir/$targetdir/$dir/Freesurfer_result/SUMA"
	echo "Current folder: $PWD"

	echo "3dcopy -prefix brainAFNI+orig brain.nii.gz"
	#3dAFNItoNIFTI -prefix brainAFNI+orig brain.nii.gz
	rm brainAFNI+orig*
	3dcopy brain.nii.gz brainAFNI+orig

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


