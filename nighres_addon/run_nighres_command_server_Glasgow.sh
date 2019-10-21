#!/bin/bash

INSTROUT=$1

if [ -z "$1" ]
then
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'Examples call from shell: (use double quotes and single quotes like in the example). For the brain.mgdm_segmentation volumes must be in RAI coordinate system, use 3dresample with the -orient RAI option to achive the right orientation of the volumes'
echo
echo 'run_nighres_command.sh '\''brain.mgdm_segmentation(contrast_image1="T1_orig_ss_nudge_RAI.nii.gz", contrast_type1="Mp2rage7T", contrast_image2="T1_ss_nudge_RAI.nii.gz", contrast_type2="T1map7T", save_data=True, file_name="sub001_test02", output_dir="resultsTestAdjustPriors_testSh/", adjust_intensity_priors=False, n_steps=5)'\'''
echo
echo 'run_nighres_command.sh '\''brain.extract_brain_region(segmentation="sub001_test02_mgdm_seg.nii.gz", levelset_boundary="sub001_test02_mgdm_dist.nii.gz", maximum_membership="sub001_test02_mgdm_mems.nii.gz", maximum_label="sub001_test02_mgdm_lbls.nii.gz", extracted_region="right_cerebrum", atlas_file=None, normalize_probabilities=False, estimate_tissue_densities=False, partial_volume_distance=1.0, save_data=True, output_dir="right_cerebrum/", file_name="sub001_right_cer")'\'''
echo
echo 'run_nighres_command.sh '\''cortex.cruise_cortex_extraction(init_image="sub001_right_cer_xmask_rcrwm.nii.gz", wm_image="sub001_right_cer_xproba_rcrwm.nii.gz", gm_image="sub001_right_cer_xproba_rcrgm.nii.gz", csf_image="sub001_right_cer_xproba_rcrbg.nii.gz", vd_image=None, data_weight=0.4, regularization_weight=0.1, max_iterations=500, normalize_probabilities=False, correct_wm_pv=True, wm_dropoff_dist=1.0, topology="wcs", topology_lut_dir=None, save_data=True, output_dir="cruise_right/", file_name="subj001_cruise_right")'\'''
echo
echo 'run_nighres_command.sh '\''laminar.volumetric_layering(inner_levelset="", outer_levelset="", n_layers=4, topology_lut_dir=None, save_data=False, output_dir="", file_name="")'\'''
echo
echo 'run_nighres_command.sh '\''laminar.profile_sampling(profile_surface_image="", intensity_image="", save_data=False, output_dir="", file_name="")'\'''
echos
echo 'run_nighres_command.sh '\''surface.probability_to_levelset(probability_image="", save_data=False, output_dir="", file_name="")'\'''
echo
echo 'run_nighres_command.sh '\''filtering.filter_ridge_structures(input_image="", structure_intensity="bright", output_type="probability", use_strict_min_max_filter=True, save_data=False, output_dir="", file_name="")'\'''
echo
echo 'run_nighres_command.sh '\''filtering.filter_ridge_structures(input_image="", structure_intensity="dark", output_type="probability", use_strict_min_max_filter=True, save_data=False, output_dir="", file_name="")'\'''
echo
exit 1
fi

INSTRCOMP=$(printf 'import sys; sys.path.insert(0,'\''%s'\''); print(sys.path); import nighres; nighres.%s' "$NIGHRES_TOOLBOXDIR" "$INSTROUT")

echo $INSTRCOMP >> _tttt.py

python _tttt.py

rm _tttt.py


