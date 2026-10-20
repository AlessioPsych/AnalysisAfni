#!/bin/bash

# Instructions to run:
# cd /media/alessiofracasso/DATADRIVE1/Flanker/code
# sh 09_00_getRoisInIndividualSpace.sh

maindir="/home/fracasso/Data/openNeuro/ds001751_FlankerContextualModulation"
codedir="code"
Freesurferdir="derivatives/Freesurfer_output"
sumaMNIDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas"
bensonRetinotopy='derivatives/surfaceAtlases/MNI_Glasser_HCP_v1.0/osfstorage-archive/prfresultsnii/999999_copy'

echo "main folder:"
echo "$maindir"
echo "code folder:"
echo "$codedir"
echo "Freesurfer folder:"
echo "$Freesurferdir"
echo "Benson Retinotopy"
echo "$bensonRetinotopy"

cd $maindir

\ls -d sub-* > subjList.txt

while IFS= read -r dir; do

    echo "Processing directory: $dir"
    echo "Freesurferdir directory: $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
        
    cd "$maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA" || { echo "Failed to enter $dir"; continue; }                
        
    echo "Current folder: $PWD"

    # Benson Retinotopy
    echo "Benson Retinotopy"
    
    # clean up
    [ -f transformed_subcortex.fit1_ang.nii.gz ] && rm transformed_subcortex.fit1_ang.nii.gz
    [ -f transformed_subcortex.fit1_ecc.nii.gz ] && rm transformed_subcortex.fit1_ecc.nii.gz
    [ -f transformed_subcortex.fit1_gain.nii.gz ] && rm transformed_subcortex.fit1_gain.nii.gz
    [ -f transformed_subcortex.fit1_meanvol.nii.gz ] && rm transformed_subcortex.fit1_meanvol.nii.gz
    [ -f transformed_subcortex.fit1_R2.nii.gz ] && rm transformed_subcortex.fit1_R2.nii.gz
    [ -f transformed_subcortex.fit1_rfsize.nii.gz ] && rm transformed_subcortex.fit1_rfsize.nii.gz
    [ -f transformed_MNI_Atlas_zp.nii.gz ] && rm transformed_MNI_Atlas_zp.nii.gz
    [ -f transformed_coll_rois.nii.gz ] && rm transformed_coll_rois.nii.gz
    [ -f transformed_BensonAtlasVolume.nii.gz ] && rm transformed_BensonAtlasVolume.nii.gz
    
 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_ang.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_ang.nii.gz
 
 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_ecc.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_ecc.nii.gz
 
 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_gain.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_gain.nii.gz
 
 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_meanvol.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_meanvol.nii.gz

 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_R2.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_R2.nii.gz

 	mri_vol2vol --targ $maindir/$bensonRetinotopy/subcortex.fit1_rfsize.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_subcortex.fit1_rfsize.nii.gz

 	mri_vol2vol --targ $maindir/$bensonRetinotopy/coll_rois.nii.gz \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv --nearest \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_coll_rois.nii.gz
    
 	mri_vol2vol --targ $maindir/$bensonRetinotopy/MNI_Atlas_zp.nii \
     	--mov $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/orig.mgz \
     	--xfm $maindir/$Freesurferdir/$dir/Freesurfer_result/mri/transforms/talairach.xfm  \
     	--inv \
     	--o $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA/transformed_MNI_Atlas_zp.nii.gz

     
     	3dTcat -prefix transformed_BensonAtlasVolume.nii.gz transformed_subcortex.fit1_ang.nii.gz transformed_subcortex.fit1_ecc.nii.gz transformed_subcortex.fit1_gain.nii.gz transformed_subcortex.fit1_meanvol.nii.gz transformed_subcortex.fit1_R2.nii.gz transformed_subcortex.fit1_rfsize.nii.gz
     
         # clean up
    [ -f transformed_subcortex.fit1_ang.nii.gz ] && rm transformed_subcortex.fit1_ang.nii.gz
    [ -f transformed_subcortex.fit1_ecc.nii.gz ] && rm transformed_subcortex.fit1_ecc.nii.gz
    [ -f transformed_subcortex.fit1_gain.nii.gz ] && rm transformed_subcortex.fit1_gain.nii.gz
    [ -f transformed_subcortex.fit1_meanvol.nii.gz ] && rm transformed_subcortex.fit1_meanvol.nii.gz
    [ -f transformed_subcortex.fit1_R2.nii.gz ] && rm transformed_subcortex.fit1_R2.nii.gz
    [ -f transformed_subcortex.fit1_rfsize.nii.gz ] && rm transformed_subcortex.fit1_rfsize.nii.gz

     echo "............."
     echo "............."
     echo "............."
     echo "............."
    
done < subjList.txt


