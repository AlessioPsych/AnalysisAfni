#!/bin/bash

# to run, remember to change the participant folder in line 12 to either 
#122017_high_res/ 
#or
#152017__high_res/
#or
#122017_resample/ 
#or
#152017_resample/

mainFolder="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/"
dataFolder="122017_resample"

cd $mainFolder
cd $dataFolder

rightBoundaries="rightVolumetricBoundaries.nii"
echo "right boundaries file: $rightBoundaries"

leftBoundaries="leftVolumetricBoundaries.nii"
echo "left boundaries file: $leftBoundaries"

if [ -d "right_surfaces_folder" ]; then
	echo "removing folder right_surfaces_folder" 
	rm -R right_surfaces_folder
fi

if [ -d "left_surfaces_folder" ]; then
	echo "removing folder left_surfaces_folder" 
	rm -R left_surfaces_folder
fi

if [ -f "boundariesThr.nii" ]; then
	echo "removing file boundariesThr.nii" 
	rm boundariesThr.nii
fi

3dcalc -a $rightBoundaries -expr '( within(a,-1000,0) )' -prefix boundariesThr.nii

generateSurfacesFromBoundaries.sh boundariesThr.nii 12 0-1-2-3-5-6 250 0.95

mv surfaces_folder/ right_surfaces_folder/

if [ -f "boundariesThr.nii" ]; then
	echo "removing file boundariesThr.nii" 
	rm boundariesThr.nii
fi

3dcalc -a $leftBoundaries -expr '( within(a,-1000,0) )' -prefix boundariesThr.nii

generateSurfacesFromBoundaries.sh boundariesThr.nii 12 0-1-2-3-5-6 250 0.95

mv surfaces_folder/ left_surfaces_folder/

if [ -f "boundariesThr.nii" ]; then
	echo "removing file boundariesThr.nii" 
	rm boundariesThr.nii
fi

