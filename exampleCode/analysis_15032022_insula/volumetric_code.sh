run_nighres_command.sh 'surface.probability_to_levelset(probability_image="white_matter_mask.nii.gz", save_data=True, output_dir="whiteMatter", file_name="whiteMatterLevelset.nii.gz")'
cp whiteMatter/*.nii.gz whiteMatterLevelset.nii.gz
 
run_nighres_command.sh 'surface.probability_to_levelset(probability_image="gray_matter_mask_out.nii.gz", save_data=True, output_dir="greyMatter", file_name="greyMatterLevelset.nii.gz")'
cp greyMatter/*.nii.gz greyMatterLevelset.nii.gz

run_nighres_command.sh 'laminar.volumetric_layering(inner_levelset="whiteMatterLevelset.nii.gz", outer_levelset="greyMatterLevelset.nii.gz", n_layers=7, topology_lut_dir=None, save_data=True, output_dir="volumetric", file_name="volumetricData.nii.gz")'

cp volumetric/volumetricData_layering-boundaries.nii.gz volumetricData_layering_boundaries.nii.gz
cp volumetric/volumetricData_layering-depth.nii.gz volumetricData_layering_depth.nii.gz
cp volumetric/volumetricData_layering-layers.nii.gz volumetricData_layering_layers.nii.gz

