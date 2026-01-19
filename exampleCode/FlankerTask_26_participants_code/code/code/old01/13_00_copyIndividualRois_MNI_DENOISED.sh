#!/bin/bash

# Instructions to run:
# cd /media/alessiofracasso/DATADRIVE1/Flanker/code
# sh 13_00_copyIndividualRois_MNI_DENOISED.sh

maindir="/media/alessiofracasso/DATADRIVE1/Flanker"
codedir="/code"
Freesurferdir="derivatives/Freesurfer_output"
inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised"
inputAfniNoDenoisedOrigDir="derivatives/processing_afni_no_denoised"
inputAfniDenoisedMNIDir="derivatives/processing_afni_MNI_denoised"
inputAfniNoDenoisedMNIDir="derivatives/processing_afni_MNI_no_denoised"
sumaMNIDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas"
glasserDir="derivatives/surfaceAtlases/MNI_Glasser_HCP_v1.0/MNI_Glasser_HCP_2019_v1.0"
princetonDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas/probatlas_v4/ProbAtlas_v4/subj_surf_all"

echo "main folder:"
echo "$maindir"
echo "code folder:"
echo "$codedir"
echo "Freesurfer folder:"
echo "$Freesurferdir"
echo "input afni denoised orig folder:"
echo "$inputAfniDenoisedOrigDir"
echo "input afni no denoised orig folder:"
echo "$inputAfniNoDenoisedOrigDir"
echo "input afni denoised MNI folder:"
echo "$inputAfniDenoisedMNIDir"
echo "input afni no denoised MNI folder:"
echo "$inputAfniNoDenoisedMNIDir"
echo "Glasser 2016 atlas folder:"
echo "$glasserDir"
echo "Wang 2015 atlas folder:"
echo "$princetonDir"

cd $maindir

\ls -d sub-* > subjList.txt

while IFS= read -r dir; do

    freeDirTemp="$sumaMNIDir"
    afniDirTempDenoised_ORIG="$maindir/$inputAfniDenoisedMNIDir/$dir.results"
    #afniDirTempNoDenoised_ORIG="$maindir/$inputAfniNoDenoisedOrigDir/$dir.results"
    #afniDirTempDenoised_MNI="$maindir/$inputAfniDenoisedMNIDir/$dir.result"
    #afniDirTempNoDenoised_MNI="$maindir/$inputAfniNoDenoisedMNIDir/$dir.results"
 
    echo "Processing directory: $dir"
    echo "Freesurferdir directory: $freeDirTemp"
    echo "afni denoised ORIG directory: $afniDirTempDenoised_ORIG" 

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz"
    echo "[ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz
        
    # copy Glasser 2016, MNI space
    echo "cp $freeDirTemp/lh.Glasser_HCP_MNI.nii.gz $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz"
    echo "cp $freeDirTemp/rh.Glasser_HCP_MNI.nii.gz $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz"
    cp $freeDirTemp/lh.Glasser_HCP_MNI.nii.gz $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz
    cp $freeDirTemp/rh.Glasser_HCP_MNI.nii.gz $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz"
    echo "[ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz

    # copy Wang 2015, MNI space
    echo "cp $freeDirTemp/lh.Wang_2015_MNI.nii.gz $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz"
    echo "cp $freeDirTemp/rh.Wang_2015_MNI.nii.gz $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz"
    cp $freeDirTemp/lh.Wang_2015_MNI.nii.gz $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz
    cp $freeDirTemp/rh.Wang_2015_MNI.nii.gz $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii"
    [ -f $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii
    
    # copy aparc+aseg, MNI space
    echo "cp $freeDirTemp/aparc+aseg.nii $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii"
    cp $freeDirTemp/aparc+aseg.nii $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii"
    [ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii
    
    # copy aparc+aseg, MIN space
    echo "cp $freeDirTemp/aparc.a2009s+aseg.nii $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii"
    cp $freeDirTemp/aparc.a2009s+aseg.nii $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii
    
    # clean up
        [ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/aparc+aseg_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_MNI_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_epiResampled_$dir.nii.gz

    # resample rois
    3dresample -input $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/lh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc
    3dresample -input $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/rh.Glasser_HCP_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc
    3dresample -input $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/lh.Wang_2015_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc
    3dresample -input $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/rh.Wang_2015_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc
    3dresample -input $afniDirTempDenoised_ORIG/aparc+aseg_MNI_$dir.nii -rmode NN -prefix $afniDirTempDenoised_ORIG/aparc+aseg_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc
    3dresample -input $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_$dir.nii -rmode NN -prefix $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_MNI_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+tlrc

    echo "............."
    echo "............."
    echo "............."
    echo "............."
    
done < subjList.txt


