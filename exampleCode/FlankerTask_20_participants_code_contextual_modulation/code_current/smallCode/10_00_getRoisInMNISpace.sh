#!/bin/bash

maindir="/media/alessiofracasso/DATADRIVE1/Flanker"
codedir="/code"
#Freesurferdir="derivatives/Freesurfer_output"
#inputAfniDenoisedOrigDir="derivatives/processing_afni_denoised"
#inputAfniNoDenoisedOrigDir="derivatives/processing_afni_no_denoised"
#inputAfniDenoisedMNIDir="derivatives/processing_afni_MNI_denoised"
#inputAfniNoDenoisedMNIDir="derivatives/processing_afni_MNI_no_denoised"
sumaMNIDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas"
glasserDir="derivatives/surfaceAtlases/MNI_Glasser_HCP_v1.0/MNI_Glasser_HCP_2019_v1.0"
princetonDir="derivatives/surfaceAtlases/suma_MNI152_2009_princetonAtlas/probatlas_v4/ProbAtlas_v4/subj_surf_all"

echo "main folder:"
echo "$maindir"
echo "code folder:"
echo "$codedir"
#echo "Freesurfer folder:"
#echo "$Freesurferdir"
#echo "input afni denoised orig folder:"
#echo "$inputAfniDenoisedOrigDir"
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

cd $maindir
cd $sumaMNIDir

echo "Current folder: $PWD"

# Glasser Atlas to MNI volume
echo "Glasser (2016) atlas to MNI volume...."

# clean up
[ -f lh.std.141.Glasser_HCP.lbl.niml.dset ] && rm lh.std.141.Glasser_HCP.lbl.niml.dset
[ -f lh.std.141.Glasser_HCP.niml.dset ] && rm lh.std.141.Glasser_HCP.niml.dset
[ -f rh.std.141.Glasser_HCP.lbl.niml.dset ] && rm rh.std.141.Glasser_HCP.lbl.niml.dset
[ -f rh.std.141.Glasser_HCP.niml.dset ] && rm rh.std.141.Glasser_HCP.niml.dset

# copy surface files Glasser 2016    
echo "cp $maindir/$glasserDir/lh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$sumaMNIDir"
echo "cp $maindir/$glasserDir/lh.std.141.Glasser_HCP.niml.dset $maindir/$sumaMNIDir"
echo "cp $maindir/$glasserDir/rh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$sumaMNIDir"
echo "cp $maindir/$glasserDir/rh.std.141.Glasser_HCP.niml.dset $maindir/$sumaMNIDir"
cp $maindir/$glasserDir/lh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$sumaMNIDir
cp $maindir/$glasserDir/lh.std.141.Glasser_HCP.niml.dset $maindir/$sumaMNIDir
cp $maindir/$glasserDir/rh.std.141.Glasser_HCP.lbl.niml.dset $maindir/$sumaMNIDir
cp $maindir/$glasserDir/rh.std.141.Glasser_HCP.niml.dset $maindir/$sumaMNIDir

# clean up
[ -f lh.Glasser_HCP_MNI.nii.gz ] && rm lh.Glasser_HCP_MNI.nii.gz
[ -f rh.Glasser_HCP_MNI.nii.gz ] && rm rh.Glasser_HCP_MNI.nii.gz

@surf_to_vol_spackle -spec std.141.MNI152_2009_lh.spec -surfA std.141.lh.smoothwm.gii -surfB std.141.lh.pial.gii -surfset lh.std.141.Glasser_HCP.niml.dset -mode -prefix lh.Glasser_HCP_MNI -maskset MNI152_2009_SurfVol.nii -meanrad 0.7 -maxiters 1

@surf_to_vol_spackle -spec std.141.MNI152_2009_rh.spec -surfA std.141.rh.smoothwm.gii -surfB std.141.rh.pial.gii -surfset rh.std.141.Glasser_HCP.niml.dset -mode -prefix rh.Glasser_HCP_MNI -maskset MNI152_2009_SurfVol.nii -meanrad 0.7 -maxiters 1

# Wang (2015) Atlas to MNI volume
echo "Wang (2015) atlas to MNI volume...."
    
# clean up
[ -f maxprob_surf_lh.1D.dset ] && rm maxprob_surf_lh.1D.dset
[ -f maxprob_surf_rh.1D.dset ] && rm maxprob_surf_rh.1D.dset

# copy surface files Wang 2015    
echo "cp $maindir/$princetonDir/maxprob_surf_lh.1D.dset $maindir/$sumaMNIDir"
echo "cp $maindir/$princetonDir/maxprob_surf_rh.1D.dset $maindir/$sumaMNIDir"
cp $maindir/$princetonDir/maxprob_surf_lh.1D.dset $maindir/$sumaMNIDir
cp $maindir/$princetonDir/maxprob_surf_rh.1D.dset $maindir/$sumaMNIDir

# clean up
[ -f lh.Wang_2015_MNI.nii.gz ] && rm lh.Wang_2015_MNI.nii.gz
[ -f rh.Wang_2015_MNI.nii.gz ] && rm rh.Wang_2015_MNI.nii.gz

@surf_to_vol_spackle -spec std.141.MNI152_2009_lh.spec -surfA std.141.lh.smoothwm.gii -surfB std.141.lh.pial.gii -surfset maxprob_surf_lh.1D.dset -mode -prefix lh.Wang_2015_MNI -maskset MNI152_2009_SurfVol.nii -meanrad 0.7 -maxiters 1

@surf_to_vol_spackle -spec std.141.MNI152_2009_rh.spec -surfA std.141.rh.smoothwm.gii -surfB std.141.rh.pial.gii -surfset maxprob_surf_rh.1D.dset -mode -prefix rh.Wang_2015_MNI -maskset MNI152_2009_SurfVol.nii -meanrad 0.7 -maxiters 1



