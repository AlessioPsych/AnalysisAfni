#!/bin/bash

# Instructions to run:
# sh 01_01_copyAnatInDenoisedFolder.sh

maindir="/media/alessiofracasso/DATADRIVE1/Flanker"
targetdir="derivatives/mrtrix3"

echo "$maindir"
echo "$targetdir"

cd "$maindir"
echo "Current folder: $PWD"

# print list of participants
\ls -d sub-* > subdirs.txt

while IFS= read -r dir; do

	cd "$maindir"

	echo "Processing directory: $dir"
	
	if [ -d $maindir/$targetdir/$dir/anat ]; then
	        echo "removing folder $maindir/$targetdir/$dir/anat" 
		rm -R $maindir/$targetdir/$dir/anat
	fi
	
	echo "data folder: $maindir/$dir/anat"
	echo "target folder: $maindir/$targetdir/$dir/anat"
	
	# copy anat folder to denoise
	echo "cp -R $maindir/$dir/anat $maindir/$targetdir/$dir/anat"	
	cp -R $maindir/$dir/anat $maindir/$targetdir/$dir/anat
	

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


