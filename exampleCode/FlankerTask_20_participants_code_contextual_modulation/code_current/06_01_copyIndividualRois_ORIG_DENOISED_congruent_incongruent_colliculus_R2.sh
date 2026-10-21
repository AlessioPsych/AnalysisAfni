#!/bin/bash

# Instructions to run:
# cd /media/alessiofracasso/DATADRIVE1/Flanker/code
# sh 11_00_copyIndividualRoisInOrigSpace.sh

maindir="/home/fracasso/Data/openNeuro/ds001751_FlankerContextualModulation"
codedir="code"
Freesurferdir="derivatives/Freesurfer_output"
inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised_incongruent_congruent_tent"
sumaMNIDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas"
glasserDir="derivatives/surfaceAtlases/MNI_Glasser_HCP_v1.0/MNI_Glasser_HCP_2019_v1.0"
princetonDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas/probatlas_v4/ProbAtlas_v4/subj_surf_all"
bensonRetinotopy='derivatives/surfaceAtlases/MNI_Glasser_HCP_v1.0/inHouse/participants_fit_Benson_jov/SUMA'

echo "main folder:"
echo "$maindir"
echo "code folder:"
echo "$codedir"
echo "Freesurfer folder:"
echo "$Freesurferdir"
echo "input afni denoised orig folder:"
echo "$inputAfniDenoisedOrigDir"
echo "Glasser 2016 atlas folder:"
echo "$glasserDir"
echo "Wang 2015 atlas folder:"
echo "$princetonDir"
echo "Benson Retinotopy folder:"
echo "$bensonRetinotopy"

cd $maindir

\ls -d sub-* > subjList.txt

while IFS= read -r dir; do

    freeDirTemp="$maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    afniDirTempDenoised_ORIG="$maindir/$inputAfniDenoisedOrigDir/$dir.results"
    #afniDirTempNoDenoised_ORIG="$maindir/$inputAfniNoDenoisedOrigDir/$dir.results"
    #afniDirTempDenoised_MNI="$maindir/$inputAfniDenoisedMNIDir/$dir.result"
    #afniDirTempNoDenoised_MNI="$maindir/$inputAfniNoDenoisedMNIDir/$dir.results"
 
    echo "Processing directory: $dir"
    echo "Freesurferdir directory: $freeDirTemp"
    echo "afni denoised ORIG directory: $afniDirTempDenoised_ORIG" 
  
    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz
        
    # copy coll ROIs, ORIG space
    echo "cp $freeDirTemp/transformed_coll_rois.nii.gz $afniDirTempDenoised_ORIG/transformed_coll_rois.nii.gz"
    cp $freeDirTemp/transformed_coll_rois.nii.gz $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz


    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_$dir.nii.gz"
	[ -f $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG transformed_BensonAtlasVolume_$dir.nii.gz        
    
    # copy Benson model results, ORIG space
    echo "cp $freeDirTemp/transformed_BensonAtlasVolume.nii.gz $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume.nii.gz"
    cp $freeDirTemp/transformed_BensonAtlasVolume.nii.gz $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_$dir.nii.gz
    
    # clean up
    [ -f $afniDirTempDenoised_ORIG/transformed_coll_rois_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/transformed_coll_rois_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_epiResampled_$dir.nii.gz
    

    # resample rois
    3dresample -input $afniDirTempDenoised_ORIG/transformed_coll_rois_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/transformed_coll_rois_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_$dir.nii.gz -rmode Lin -prefix $afniDirTempDenoised_ORIG/transformed_BensonAtlasVolume_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig

    echo "............."
    echo "............."
    echo "............."
    echo "............."
    
done < subjList.txt


