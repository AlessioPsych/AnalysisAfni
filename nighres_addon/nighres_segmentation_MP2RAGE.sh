#!/bin/bash

T1MP2RAGE_SS_RAI=$1
T1MAP_SS_RAI=$2

if [ -z "$1" ]
then
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'T1MP2RAGE_SS_RAI=$1, MP2RAGE T1 contrast, skull stripped in RAI coordinates (check 3dresample with -orient RAI flag)'
echo 'T1MAP_SS_RAI=$2, T1 map, skull stripped, in RAI coordinates' 
echo
exit 1
fi

cp $T1MP2RAGE_SS_RAI _ttt_T1_MP2RAGE_RAI.nii.gz
cp $T1MAP_SS_RAI _ttt_T1_MAP_RAI.nii.gz

run_nighres_command.sh 'brain.mgdm_segmentation(contrast_image1="_ttt_T1_MP2RAGE_RAI.nii.gz", contrast_type1="T1map7T", save_data=True, file_name="MP2RAGE_seg", output_dir="nighres_resultsMP2RAGE_T1M/", n_steps=5, max_iterations=800)'


#run_nighres_command.sh 'brain.mgdm_segmentation(contrast_image1="_ttt_T1_MP2RAGE_RAI.nii.gz", contrast_type1="Mp2rage7T", save_data=True, #file_name="MP2RAGE_seg", output_dir="nighres_resultsMP2RAGE/", n_steps=5, max_iterations=800, atlas_file="/home/fracasso/programs/nighres/nighres/atlases/#brain-segmentation-prior3.0/brain-atlas-3.0.3_mod.txt")'

#run_nighres_command.sh 'brain.mgdm_segmentation(contrast_image1="_ttt_T1_MP2RAGE_RAI.nii.gz", contrast_type1="Mp2rage7T", contrast_image2="_ttt_T1_MAP_RAI.nii.gz", contrast_type2="T1map7T", save_data=True, file_name="MP2RAGE_seg", output_dir="nighres_resultsMP2RAGE/", n_steps=5, max_iterations=800, atlas_file="/home/fracasso/programs/nighres/nighres/atlases/brain-segmentation-prior3.0/brain-atlas-3.0.3_mod.txt")'


cd nighres_resultsMP2RAGE_T1M

run_nighres_command.sh 'brain.extract_brain_region(segmentation="MP2RAGE_seg_mgdm_seg.nii.gz", levelset_boundary="MP2RAGE_seg_mgdm_dist.nii.gz", maximum_membership="MP2RAGE_seg_mgdm_mems.nii.gz", maximum_label="MP2RAGE_seg_mgdm_lbls.nii.gz", extracted_region="right_cerebrum", atlas_file=None, normalize_probabilities=False, estimate_tissue_densities=False, partial_volume_distance=1.0, save_data=True, output_dir="right_cerebrum/", file_name="MP2RAGE_seg_right_cer")'

run_nighres_command.sh 'brain.extract_brain_region(segmentation="MP2RAGE_seg_mgdm_seg.nii.gz", levelset_boundary="MP2RAGE_seg_mgdm_dist.nii.gz", maximum_membership="MP2RAGE_seg_mgdm_mems.nii.gz", maximum_label="MP2RAGE_seg_mgdm_lbls.nii.gz", extracted_region="left_cerebrum", atlas_file=None, normalize_probabilities=False, estimate_tissue_densities=False, partial_volume_distance=1.0, save_data=True, output_dir="left_cerebrum/", file_name="MP2RAGE_seg_left_cer")'

cd right_cerebrum

run_nighres_command.sh 'cortex.cruise_cortex_extraction(init_image="MP2RAGE_seg_right_cer_xmask_rcrwm.nii.gz", wm_image="MP2RAGE_seg_right_cer_xproba_rcrwm.nii.gz", gm_image="MP2RAGE_seg_right_cer_xproba_rcrgm.nii.gz", csf_image="MP2RAGE_seg_right_cer_xproba_rcrbg.nii.gz", vd_image=None, data_weight=0.4, regularization_weight=0.1, max_iterations=500, normalize_probabilities=False, correct_wm_pv=True, wm_dropoff_dist=1.0, topology="wcs", topology_lut_dir=None, save_data=True, output_dir="cruise_right/", file_name="MP2RAGE_seg_cruise_right")'

cd ..
cd left_cerebrum

run_nighres_command.sh 'cortex.cruise_cortex_extraction(init_image="MP2RAGE_seg_left_cer_xmask_lcrwm.nii.gz", wm_image="MP2RAGE_seg_left_cer_xproba_lcrwm.nii.gz", gm_image="MP2RAGE_seg_left_cer_xproba_lcrgm.nii.gz", csf_image="MP2RAGE_seg_left_cer_xproba_lcrbg.nii.gz", vd_image=None, data_weight=0.4, regularization_weight=0.1, max_iterations=500, normalize_probabilities=False, correct_wm_pv=True, wm_dropoff_dist=1.0, topology="wcs", topology_lut_dir=None, save_data=True, output_dir="cruise_left/", file_name="MP2RAGE_seg_cruise_left")'

cd ../..

cp nighres_resultsMP2RAGE_T1M/left_cerebrum/cruise_left/MP2RAGE_seg_cruise_left_cruise_cortex.nii.gz nighres_resultsMP2RAGE_T1M/ 
cp nighres_resultsMP2RAGE_T1M/right_cerebrum/cruise_right/MP2RAGE_seg_cruise_right_cruise_cortex.nii.gz nighres_resultsMP2RAGE_T1M/ 
cp nighres_resultsMP2RAGE_T1M/left_cerebrum/cruise_left/MP2RAGE_seg_cruise_left_cruise_thick.nii.gz nighres_resultsMP2RAGE_T1M/ 
cp nighres_resultsMP2RAGE_T1M/right_cerebrum/cruise_right/MP2RAGE_seg_cruise_right_cruise_thick.nii.gz nighres_resultsMP2RAGE_T1M/ 

rm _ttt_T1_MP2RAGE_RAI.nii.gz
rm _ttt_T1_MAP_RAI.nii.gz
