
run_nighres_command.sh 'laminar.profile_sampling(profile_surface_image="volumetricData_layering_boundaries.nii.gz", intensity_image="inputData.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii.gz")'

cp profilesDir/profilesOut** profilesOut_profiles.nii.gz
rm -R profilesDir
