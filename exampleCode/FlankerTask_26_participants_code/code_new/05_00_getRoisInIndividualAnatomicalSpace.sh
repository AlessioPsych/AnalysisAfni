#!/bin/bash

#to run ' sh 05_00_getRoisInIndividualAnatomicalSpace.sh '


maindir="/mnt/disk01/ds000102_FlankerTask"
codedir="code_new"
Freesurferdir="derivatives/Freesurfer_output"
inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised_incongruent_congruent"
#inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised"
#inputAfniNoDenoisedOrigDir="derivatives/processing_afni_no_denoised"
#inputAfniDenoisedMNIDir="derivatives/processing_afni_MNI_denoised"
#inputAfniNoDenoisedMNIDir="derivatives/processing_afni_MNI_no_denoised"
sumaMNIDir="/mnt/disk01/surfaceAtlas_HCP_vonEconomo/surfaceAtlases/suma_MNI152_2009_princetonAtlas"
glasserDir="/mnt/disk01/surfaceAtlas_HCP_vonEconomo/surfaceAtlases/MNI_Glasser_HCP_v1.0/MNI_Glasser_HCP_2019_v1.0"
princetonDir="/mnt/disk01/surfaceAtlas_HCP_vonEconomo/surfaceAtlases/suma_MNI152_2009_princetonAtlas/probatlas_v4/ProbAtlas_v4/subj_surf_all"
bensonRetinotopy='/mnt/disk01/surfaceAtlas_HCP_vonEconomo/surfaceAtlases/MNI_Glasser_HCP_v1.0/inHouse/participants_fit_Benson_jov/SUMA'

echo "main folder:"
echo "$maindir"
echo "code folder:"
echo "$codedir"
echo "Freesurfer folder:"
echo "$Freesurferdir"
echo "input afni denoised orig folder:"
echo "$inputAfniDenoisedOrigDir"
#echo "input afni no denoised orig folder:"
#echo "$inputAfniNoDenoisedOrigDir"
#echo "input afni denoised MNI folder:"
#echo "$inputAfniDenoisedMNIDir"
#echo "input afni no denoised MNI folder:"
#echo "$inputAfniNoDenoisedMNIDir"
echo "Glasser 2016 atlas folder:"
echo "$glasserDir"
echo "Wang 2015 atlas folder:"
echo "$princetonDir"
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
    [ -f std.141.lh.999999.1D.dset ] && rm std.141.lh.999999.1D.dset
    [ -f std.141.rh.999999.1D.dset ] && rm std.141.rh.999999.1D.dset
    
    echo "cp $bensonRetinotopy/std.141.lh.999999.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    echo "cp $bensonRetinotopy/std.141.rh.999999.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    cp $bensonRetinotopy/std.141.lh.999999.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA
    cp $bensonRetinotopy/std.141.rh.999999.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA

    # clean up
    [ -f lh.Benson_Retinotopy.nii.gz ] && rm lh.Benson_Retinotopy.nii.gz
    [ -f rh.Benson_Retinotopy.nii.gz ] && rm rh.Benson_Retinotopy.nii.gz

    @surf_to_vol_spackle -spec std.141.Freesurfer_result_lh.spec -surfA std.141.lh.smoothwm.gii -surfB std.141.lh.pial.gii -surfset std.141.lh.999999.1D.dset -prefix lh.Benson_Retinotopy -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1
    
    @surf_to_vol_spackle -spec std.141.Freesurfer_result_rh.spec -surfA std.141.rh.smoothwm.gii -surfB std.141.rh.pial.gii -surfset std.141.rh.999999.1D.dset -prefix rh.Benson_Retinotopy -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1

    # Glasser Atlas
    echo "Glasser (2016) atlas...."

    # clean up
    [ -f lh.std.141.Glasser_HCP.lbl.niml.dset ] && rm lh.std.141.Glasser_HCP.lbl.niml.dset
    [ -f lh.std.141.Glasser_HCP.niml.dset ] && rm lh.std.141.Glasser_HCP.niml.dset
    [ -f rh.std.141.Glasser_HCP.lbl.niml.dset ] && rm rh.std.141.Glasser_HCP.lbl.niml.dset
    [ -f rh.std.141.Glasser_HCP.niml.dset ] && rm rh.std.141.Glasser_HCP.niml.dset

    # copy surface files Glasser 2016    
    echo "cp $glasserDir/lh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    echo "cp $glasserDir/lh.std.141.Glasser_HCP.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    echo "cp $glasserDir/rh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    echo "cp $glasserDir/rh.std.141.Glasser_HCP.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    cp $glasserDir/lh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA
    cp $glasserDir/lh.std.141.Glasser_HCP.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA
    cp $glasserDir/rh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA
    cp $glasserDir/rh.std.141.Glasser_HCP.niml.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA

    # clean up
    [ -f lh.Glasser_HCP.nii.gz ] && rm lh.Glasser_HCP.nii.gz
    [ -f rh.Glasser_HCP.nii.gz ] && rm rh.Glasser_HCP.nii.gz

    @surf_to_vol_spackle -spec std.141.Freesurfer_result_lh.spec -surfA std.141.lh.smoothwm.gii -surfB std.141.lh.pial.gii -surfset lh.std.141.Glasser_HCP.niml.dset -mode -prefix lh.Glasser_HCP -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1

    @surf_to_vol_spackle -spec std.141.Freesurfer_result_rh.spec -surfA std.141.rh.smoothwm.gii -surfB std.141.rh.pial.gii -surfset rh.std.141.Glasser_HCP.niml.dset -mode -prefix rh.Glasser_HCP -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1
          
    # Wang (2015) Atlas
    echo "Wang (2015) atlas...."
    
    # clean up
    [ -f maxprob_surf_lh.1D.dset ] && rm maxprob_surf_lh.1D.dset
    [ -f maxprob_surf_rh.1D.dset ] && rm maxprob_surf_rh.1D.dset

    # copy surface files Wang 2015    
    echo "cp $princetonDir/maxprob_surf_lh.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    echo "cp $princetonDir/maxprob_surf_rh.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA"
    cp $princetonDir/maxprob_surf_lh.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA
    cp $princetonDir/maxprob_surf_rh.1D.dset $maindir/$Freesurferdir/$dir/Freesurfer_result/SUMA

    # clean up
    [ -f lh.Wang_2015.nii.gz ] && rm lh.Wang_2015.nii.gz
    [ -f rh.Wang_2015.nii.gz ] && rm rh.Wang_2015.nii.gz

    @surf_to_vol_spackle -spec std.141.Freesurfer_result_lh.spec -surfA std.141.lh.smoothwm.gii -surfB std.141.lh.pial.gii -surfset maxprob_surf_lh.1D.dset -mode -prefix lh.Wang_2015 -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1

    @surf_to_vol_spackle -spec std.141.Freesurfer_result_rh.spec -surfA std.141.rh.smoothwm.gii -surfB std.141.rh.pial.gii -surfset maxprob_surf_rh.1D.dset -mode -prefix rh.Wang_2015 -maskset Freesurfer_result_SurfVol.nii -meanrad 0.7 -maxiters 1
     
     
     echo "............."
     echo "............."
     echo "............."
     echo "............."
    
done < subjList.txt


