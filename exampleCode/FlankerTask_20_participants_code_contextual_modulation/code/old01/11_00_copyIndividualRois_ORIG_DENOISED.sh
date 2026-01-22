#!/bin/bash

# Instructions to run:
# cd /media/alessiofracasso/DATADRIVE1/Flanker/code
# sh 11_00_copyIndividualRoisInOrigSpace.sh

maindir="/home/fracasso/Data/openNeuro/ds000102"
codedir="/code"
Freesurferdir="derivatives/Freesurfer_output"
inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised"
inputAfniNoDenoisedOrigDir="derivatives/processing_afni_no_denoised"
inputAfniDenoisedMNIDir="derivatives/processing_afni_MNI_denoised"
inputAfniNoDenoisedMNIDir="derivatives/processing_afni_MNI_no_denoised"
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
    echo "afni no denoised ORIG directory: $afniDirTempNoDenoised_ORIG" 
    echo "afni denoised MNI directory: $afniDirTempDenoised_MNI" 
    echo "afni no denoised MNI directory: $afniDirTempNoDenoised_MNI" 

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz"
    echo "[ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz
        
    # copy Glasser 2016, ORIG space
    echo "cp $freeDirTemp/lh.Glasser_HCP.nii.gz $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz"
    echo "cp $freeDirTemp/rh.Glasser_HCP.nii.gz $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz"
    cp $freeDirTemp/lh.Glasser_HCP.nii.gz $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz
    cp $freeDirTemp/rh.Glasser_HCP.nii.gz $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz"
    echo "[ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz

    # copy Wang 2015, ORIG space
    echo "cp $freeDirTemp/lh.Wang_2015.nii.gz $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz"
    echo "cp $freeDirTemp/rh.Wang_2015.nii.gz $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz"
    cp $freeDirTemp/lh.Wang_2015.nii.gz $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz
    cp $freeDirTemp/rh.Wang_2015.nii.gz $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz
    
    # copy aparc+aseg, ORIG space
    echo "cp $freeDirTemp/aparc+aseg.nii.gz $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz"
    cp $freeDirTemp/aparc+aseg.nii.gz $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz

    # clean up
    echo "[ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz
    
    # copy aparc+aseg, ORIG space
    echo "cp $freeDirTemp/aparc.a2009s+aseg.nii.gz $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz"
    cp $freeDirTemp/aparc.a2009s+aseg.nii.gz $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz

    # clean up Benson Retinotopy Left
        echo "[ -f $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz    
    
    # clean up Benson Retinotopy Right
    echo "[ -f $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz"
    [ -f $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz    
    
    # copy Benson retinotopy Left, ORIG space
    echo "cp $freeDirTemp/lh.Benson_Retinotopy.nii.gz $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz"
    cp $freeDirTemp/lh.Benson_Retinotopy.nii.gz $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz

    # copy Benson retinotopy Right, ORIG space
    echo "cp $freeDirTemp/rh.Benson_Retinotopy.nii.gz $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz"
    cp $freeDirTemp/rh.Benson_Retinotopy.nii.gz $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz
    
    # clean up
        [ -f $afniDirTempDenoised_ORIG/lh.Glasser_HCP_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Glasser_HCP_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Glasser_HCP_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Glasser_HCP_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/lh.Wang_2015_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Wang_2015_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Wang_2015_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Wang_2015_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/aparc+aseg_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc+aseg_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_epiResampled_$dir.nii.gz
    [ -f $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_epiResampled_$dir.nii.gz ] && rm $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_epiResampled_$dir.nii.gz
    

    # resample rois
    3dresample -input $afniDirTempDenoised_ORIG/lh.Glasser_HCP_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/lh.Glasser_HCP_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/rh.Glasser_HCP_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/rh.Glasser_HCP_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/lh.Wang_2015_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/lh.Wang_2015_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/rh.Wang_2015_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/rh.Wang_2015_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/aparc+aseg_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/aparc+aseg_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/aparc.a2009s+aseg_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/lh.Benson_Retinotopy_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig
    3dresample -input $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_$dir.nii.gz -rmode NN -prefix $afniDirTempDenoised_ORIG/rh.Benson_Retinotopy_epiResampled_$dir.nii.gz  -master $afniDirTempDenoised_ORIG/stats.$dir+orig

    echo "............."
    echo "............."
    echo "............."
    echo "............."
    
done < subjList.txt


