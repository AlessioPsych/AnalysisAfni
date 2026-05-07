##!/usr/bin/env bash

MPRAGE=$1

cp $MPRAGE anatCopy.nii.gz

3dcalc -a anatCopy.nii.gz -expr 'step(a)' -prefix _ttt_anatMask.nii.gz
3dSeg    -anat anatCopy.nii.gz    -mask _ttt_anatMask.nii.gz \
                    -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' \
                    -bias_fwhm 25 -mixfrac UNI -main_N 5 \
                    -blur_meth BFT
cp Segsy/Posterior+orig.BRIK Posterior+orig.BRIK
cp Segsy/Posterior+orig.HEAD Posterior+orig.HEAD
3dAFNItoNIFTI -prefix probMaps.nii.gz Posterior+orig
3dTcat -prefix csfProb.nii.gz 'probMaps.nii.gz[0]'
3dTcat -prefix gmProb.nii.gz 'probMaps.nii.gz[1]'
3dTcat -prefix wmProb.nii.gz 'probMaps.nii.gz[2]'
3dcalc -a wmProb.nii.gz -expr 'within(a,0.5,1)' -prefix wmInit.nii.gz

run_nighres_command.sh 'cortex.cruise_cortex_extraction(init_image="wmInit.nii.gz", wm_image="wmProb.nii.gz", gm_image="gmProb.nii.gz", csf_image="csfProb.nii.gz", vd_image=None, data_weight=0.4, regularization_weight=0.1, max_iterations=500, normalize_probabilities=False, correct_wm_pv=True, wm_dropoff_dist=1.0, topology="wcs", topology_lut_dir=None, save_data=True, output_dir="cruise_whole/", file_name="test_cruise_whole")'

cp cruise_whole/test_cruise_whole_cruise_cortex.nii.gz cruise_cortex.nii.gz

3dcalc -a cruise_cortex.nii.gz -expr 'step(a)' -prefix greyMatterCruise.nii.gz
3dcalc -a cruise_cortex.nii.gz -expr 'within(a,1.9,2.1)' -prefix whiteMatterCruise.nii.gz

run_nighres_command.sh 'surface.probability_to_levelset(probability_image="whiteMatterCruise.nii.gz", save_data=True, output_dir="whiteMatter", file_name="whiteMatterLevelset.nii.gz")'
cp whiteMatter/*.nii.gz whiteMatterLevelset.nii.gz
 
run_nighres_command.sh 'surface.probability_to_levelset(probability_image="greyMatterCruise.nii.gz", save_data=True, output_dir="greyMatter", file_name="greyMatterLevelset.nii.gz")'
cp greyMatter/*.nii.gz greyMatterLevelset.nii.gz

run_nighres_command.sh 'laminar.volumetric_layering(inner_levelset="whiteMatterLevelset.nii.gz", outer_levelset="greyMatterLevelset.nii.gz", n_layers=7, topology_lut_dir=None, save_data=True, output_dir="volumetric", file_name="volumetricData")'

cp volumetric/volumetricData_layering_boundaries.nii.gz volumetricData_layering_boundaries.nii.gz
cp volumetric/volumetricData_layering_depth.nii.gz volumetricData_layering_depth.nii.gz
cp volumetric/volumetricData_layering_layers.nii.gz volumetricData_layering_layers.nii.gz

run_nighres_command.sh 'laminar.profile_sampling(profile_surface_image="volumetricData_layering_boundaries.nii.gz", intensity_image="anatCopy.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut")'

cp profilesDir/profilesOut_profiles.nii.gz profilesOut_profiles.nii.gz



