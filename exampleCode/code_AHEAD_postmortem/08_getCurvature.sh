#!/bin/bash

# to run, remember to change the participant folder in line 12 to either 
#122017_high_res/ 
#or
#152017_high_res/

mainFolder="/media/alessiofracasso/DATADRIVE1/AHEAD_exvivo/"
dataFolder="152017_resample"

cd $mainFolder
cd $dataFolder

if [ -d "left_surfaces_folder_conv" ]; then
	echo "removing folder left_surfaces_folder_conv" 
	rm -R left_surfaces_folder_conv
fi

if [ -f "leftConvexity.nii.gz" ]; then
	echo "removing file leftConvexity.nii.gz" 
	rm leftConvexity.nii.gz
fi

if [ -d "right_surfaces_folder_conv" ]; then
	echo "removing folder right_surfaces_folder_conv" 
	rm -R right_surfaces_folder_conv
fi

if [ -f "rightConvexity.nii.gz" ]; then
	echo "removing file rightConvexity.nii.gz" 
	rm rightConvexity.nii.gz
fi

echo "compute convexity, left hemisphere"
metricSurfaceMap.sh -conv left_surfaces_folder

echo "compute convexity, right hemisphere"
metricSurfaceMap.sh -conv right_surfaces_folder

echo "project convexity on volume, left hemisphere"
3dSurf2Vol -spec left_surfaces_folder/spec.surfaces.smoothed 	\
	-surf_A left_surfaces_folder/boundary02_sm.1D.coord 	\
 	-sv PD_resampled.nii.gz					\
 	-grid_parent PD_resampled.nii.gz			\
 	-map_func ave						\
 	-prefix leftConvexity.nii.gz				\
 	-sdata_1D left_surfaces_folder_conv/boundary02_sm.metric.conv.1D.dset
 
echo "project convexity on volume, right hemisphere"
3dSurf2Vol -spec right_surfaces_folder/spec.surfaces.smoothed 	\
	-surf_A right_surfaces_folder/boundary02_sm.1D.coord 	\
 	-sv PD_resampled.nii.gz					\
 	-grid_parent PD_resampled.nii.gz			\
 	-map_func ave						\
 	-prefix rightConvexity.nii.gz				\
 	-sdata_1D right_surfaces_folder_conv/boundary02_sm.metric.conv.1D.dset

