#!/bin/bash

MPRAGE=$1
INV2=$2

if [ -z "$1" ]
then
echo 'segment script using nighres (MP2RAGE)'
echo 'Inputs:'
echo 'MPRAGE, input anatomy'
echo 'INV2, second inversion volume'
echo 'example call: segmentScript_inhouse01042018.sh mp2rage.nii.gz inv2.nii.gz'
exit 1
fi

3dresample -dxyz 0.7 0.7 0.7 -orient RAI -prefix _del_MP2RAGE_RAI.nii.gz -inset $MPRAGE -rmode Lin
3dresample -master _del_MP2RAGE_RAI.nii.gz -prefix _del_INV2_RAI.nii.gz -inset $INV2 -rmode Lin
3dSkullStrip -input _del_INV2_RAI.nii.gz -prefix _del_SKULL_STRIP.nii.gz -shrink_fac 0.4
3dcalc -a _del_SKULL_STRIP.nii.gz -b _del_MP2RAGE_RAI.nii.gz -expr 'b*step(a)' -prefix _del_MP2RAGE_RAI_SKULL_STRIP.nii.gz 
3dAutobox -input _del_MP2RAGE_RAI_SKULL_STRIP.nii.gz -noclust -prefix MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz -npad 5 
segment_white_matter_params_no_ants.sh MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz 0.0001 0.001 5 7-7-7-7-7-7 0 0 1 
segment_gray_matter.sh MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz white_matter_mask.nii.gz 0.65-0.75-0.75-0.75

cp MP2RAGE_RAI_SKULL_STRIP_BOX.nii.gz anatCopy.nii.gz

run_nighres_command.sh 'filtering.filter_ridge_structures(input_image="anatCopy.nii.gz", structure_intensity="dark", output_type="probability", use_strict_min_max_filter=True, save_data=True, output_dir="filter_dark/", file_name="filter")'

cp filter_dark/filter* filter_rdg_dark.nii.gz

run_nighres_command.sh 'filtering.filter_ridge_structures(input_image="anatCopy.nii.gz", structure_intensity="bright", output_type="probability", use_strict_min_max_filter=True, save_data=True, output_dir="filter_bright/", file_name="filter")'

cp filter_bright/filter* filter_rdg_bright.nii.gz

3dcalc -a filter_rdg_bright.nii.gz -b filter_rdg_dark.nii.gz -expr 'max(a,b)' -prefix filter_rdg.nii.gz

3dcalc -a white_matter_mask.nii.gz -b gray_matter_mask_out.nii.gz -expr 'a+b' -prefix brainMask.nii.gz

3dSeg    -anat anatCopy.nii.gz -mask brainMask.nii.gz \
                    -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' \
                    -bias_fwhm 25 -mixfrac UNI -main_N 5
cp Segsy/Posterior+orig.BRIK Posterior+orig.BRIK
cp Segsy/Posterior+orig.HEAD Posterior+orig.HEAD
3dAFNItoNIFTI -prefix probMaps.nii.gz Posterior+orig
3dTcat -prefix csfProb.nii.gz 'probMaps.nii.gz[0]'
3dTcat -prefix gmProb.nii.gz 'probMaps.nii.gz[1]'
3dTcat -prefix wmProb.nii.gz 'probMaps.nii.gz[2]'

3dcalc -a csfProb.nii.gz -expr 'step(a-0.8)*a' -prefix csfProb_filter.nii.gz

run_nighres_command.sh 'cortex.cruise_cortex_extraction(init_image="white_matter_mask.nii.gz", wm_image="white_matter_mask.nii.gz", gm_image="gray_matter_mask_out.nii.gz", csf_image="csfProb_filter.nii.gz", vd_image="filter_rdg.nii.gz", data_weight=0.4, regularization_weight=0.1, max_iterations=500, normalize_probabilities=False, correct_wm_pv=True, wm_dropoff_dist=1.0, topology="wcs", topology_lut_dir=None, save_data=True, output_dir="cruise_whole/", file_name="test_cruise_whole")'

cp cruise_whole/test_cruise_whole_cruise_cortex.nii.gz cruise_cortex.nii.gz

3dcalc -a cruise_cortex.nii.gz -expr 'step(a)' -prefix greyMatterCruise.nii.gz
3dcalc -a cruise_cortex.nii.gz -expr 'within(a,1.9,2.1)' -prefix whiteMatterCruise.nii.gz

run_nighres_command.sh 'surface.probability_to_levelset(probability_image="whiteMatterCruise.nii.gz", save_data=True, output_dir="whiteMatter", file_name="whiteMatterLevelset.nii.gz")'
cp whiteMatter/*.nii.gz whiteMatterLevelset.nii.gz
 
run_nighres_command.sh 'surface.probability_to_levelset(probability_image="greyMatterCruise.nii.gz", save_data=True, output_dir="greyMatter", file_name="greyMatterLevelset.nii.gz")'
cp greyMatter/*.nii.gz greyMatterLevelset.nii.gz

run_nighres_command.sh 'laminar.volumetric_layering(inner_levelset="whiteMatterLevelset.nii.gz", outer_levelset="greyMatterLevelset.nii.gz", n_layers=7, topology_lut_dir=None, save_data=True, output_dir="volumetric", file_name="volumetricData")'

cp volumetric/volumetricData_layering-boundaries.nii.gz volumetricData_layering_boundaries.nii.gz
cp volumetric/volumetricData_layering-depth.nii.gz volumetricData_layering_depth.nii.gz
cp volumetric/volumetricData_layering-layers.nii.gz volumetricData_layering_layers.nii.gz

run_nighres_command.sh 'laminar.profile_sampling(profile_surface_image="volumetricData_layering_boundaries.nii.gz", intensity_image="anatCopy.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut")'

cp profilesDir/profilesOut** profilesOut_profiles.nii.gz

#hemisphereSeparation.sh anatCopy.nii.gz cruise_cortex.nii.gz /usr/local/bin/TT_icbm452+tlrc 1
# copy script to generate two brain surfaces, in NL server;




