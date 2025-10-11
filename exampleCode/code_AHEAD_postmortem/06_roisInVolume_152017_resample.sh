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
dataFolder="152017_resample"

cd $mainFolder
cd $dataFolder

if [ -f "leftInsulaRoi.nii.gz" ]; then
	echo "removing folder leftInsulaRoi.nii.gz" 
	rm leftInsulaRoi.nii.gz
fi

if [ -f "rightInsulaRoi.nii.gz" ]; then
	echo "removing folder rightInsulaRoi.nii.gz" 
	rm rightInsulaRoi.nii.gz
fi

surf2vol.sh leftInsulaRoi.1D.roi boundary01 leftVolumetricDepth.nii 20 left_surfaces_folder/
3dAFNItoNIFTI -prefix leftInsulaRoi.nii.gz leftInsulaRoi.1D.roi_clust+orig
echo "clean up..." 
rm leftInsulaRoi.1D.roi_clust+orig*

surf2vol.sh rightInsulaRoi.1D.roi boundary02 rightVolumetricDepth.nii 20 right_surfaces_folder/
3dAFNItoNIFTI -prefix rightInsulaRoi.nii.gz rightInsulaRoi.1D.roi_clust+orig
echo "clean up..." 
rm rightInsulaRoi.1D.roi_clust+orig*

