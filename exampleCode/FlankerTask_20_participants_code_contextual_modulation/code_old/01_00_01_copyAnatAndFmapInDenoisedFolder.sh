#!/bin/bash

# Instructions to run:
# sh 01_00_01_copyAnatAndFmapInDenoisedFolder.sh

maindir="/mnt/disk01/ds001751_FlankerTask_context"
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
	
	#if [ -d $maindir/$targetdir/$dir/anat ]; then
	#        echo "removing folder $maindir/$targetdir/$dir/anat" 
	#	rm -R $maindir/$targetdir/$dir/anat
	#fi
	
	echo "data folder: $maindir/$dir/anat"
	echo "target folder: $maindir/$targetdir/$dir/anat"

	echo "data folder: $maindir/$dir/fmap"
	echo "target folder: $maindir/$targetdir/$dir/fmap"
	
	# copy anat folder to denoise
	echo "cp -R $maindir/$dir/anat $maindir/$targetdir/$dir/anat"	
	cp -R $maindir/$dir/anat $maindir/$targetdir/$dir/anat
	
	# copy fmap folder to denoise	
	echo "cp -R $maindir/$dir/fmap $maindir/$targetdir/$dir/fmap"	
	cp -R $maindir/$dir/fmap $maindir/$targetdir/$dir/fmap
	

	echo "........."
	echo "........."
	echo "........."
          
done < subdirs.txt


