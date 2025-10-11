#!/bin/bash

# to run, remember to change the participant folder in line 9, either 122017_resample/ or 152017_resample/
# sudo docker run -v /media/alessiofracasso/DATADRIVE1/AHEAD_exvivo:/dataNighres --rm -it nighres
# cd /dataNighres
# cd code/
# 

maindir="/dataNighres/"
targetdir="152017_high_res/"

cd $maindir
cd $targetdir

cruiseLeftWM="cruiseLeft_whiteMatter.nii.gz"
echo "cruise left file: $cruiseLeftWM"

cruiseLeftGM="cruiseLeft_brainMask.nii.gz"
echo "cruise left file: $cruiseLeftGM"

cruiseRightWM="cruiseRight_whiteMatter.nii.gz"
echo "cruise rught file: $cruiseRightWM"

cruiseRightGM="cruiseRight_brainMask.nii.gz"
echo "cruise right file: $cruiseRightGM"

Biel="Bieloschowsky_resampled.nii.gz"
echo "Bieloschowsky-interpolated file: $Biel"

qR1="qR1_resampled.nii.gz"
echo "quantitative R1 file: $qR1"

qR2="qR2_resampled.nii.gz"
echo "quantitative R2 file: $qR2"

PD="PD_resampled.nii.gz"
echo "PD file: $PD"

python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="cruiseLeft_whiteMatter.nii.gz", save_data=True, output_dir="cruiseLeft_whiteMatter_folder", file_name="cruiseLeft_whiteMatter_levelset.nii.gz")'
chmod 777 -R cruiseLeft_whiteMatter_folder
cp cruiseLeft_whiteMatter_folder/cruiseLeft_whiteMatter_levelset_p2l-surf.nii.gz cruiseLeft_whiteMatter_levelset_p2l-surf.nii.gz

python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="cruiseRight_whiteMatter.nii.gz", save_data=True, output_dir="cruiseRight_whiteMatter_folder", file_name="cruiseRight_whiteMatter_levelset.nii.gz")'
chmod 777 -R cruiseRight_whiteMatter_folder
cp cruiseRight_whiteMatter_folder/cruiseRight_whiteMatter_levelset_p2l-surf.nii.gz cruiseRight_whiteMatter_levelset_p2l-surf.nii.gz

python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="cruiseLeft_brainMask.nii.gz", save_data=True, output_dir="cruiseLeft_brainMask_folder", file_name="cruiseLeft_brainMask_levelset.nii.gz")'
chmod 777 -R cruiseLeft_brainMask_folder
cp cruiseLeft_brainMask_folder/cruiseLeft_brainMask_levelset_p2l-surf.nii.gz cruiseLeft_brainMask_levelset_p2l-surf.nii.gz

python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="cruiseRight_brainMask.nii.gz", save_data=True, output_dir="cruiseRight_brainMask_folder", file_name="cruiseRight_brainMask_levelset.nii.gz")'
chmod 777 -R cruiseRight_brainMask_folder
cp cruiseRight_brainMask_folder/cruiseRight_brainMask_levelset_p2l-surf.nii.gz cruiseRight_brainMask_levelset_p2l-surf.nii.gz

python3 -c 'import nighres; nighres.laminar.volumetric_layering(inner_levelset="cruiseLeft_whiteMatter_levelset_p2l-surf.nii.gz", outer_levelset="cruiseLeft_brainMask_levelset_p2l-surf.nii.gz", n_layers=7, topology_lut_dir=None, save_data=True, output_dir="leftVolumetric", file_name="leftVolumetricData.nii")'

mv leftVolumetric/leftVolumetricData_layering-boundaries.nii leftVolumetricBoundaries.nii
mv leftVolumetric/leftVolumetricData_layering-depth.nii leftVolumetricDepth.nii
mv leftVolumetric/leftVolumetricData_layering-layers.nii leftVolumetricLayers.nii

python3 -c 'import nighres; nighres.laminar.volumetric_layering(inner_levelset="cruiseRight_whiteMatter_levelset_p2l-surf.nii.gz", outer_levelset="cruiseRight_brainMask_levelset_p2l-surf.nii.gz", n_layers=7, topology_lut_dir=None, save_data=True, output_dir="rightVolumetric", file_name="rightVolumetricData.nii")'

mv rightVolumetric/rightVolumetricData_layering-boundaries.nii rightVolumetricBoundaries.nii
mv rightVolumetric/rightVolumetricData_layering-depth.nii rightVolumetricDepth.nii
mv rightVolumetric/rightVolumetricData_layering-layers.nii rightVolumetricLayers.nii

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="leftVolumetricBoundaries.nii", intensity_image="qR1_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii left_profiles_qR1.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="rightVolumetricBoundaries.nii", intensity_image="qR1_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii right_profiles_qR1.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="leftVolumetricBoundaries.nii", intensity_image="qR2_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii left_profiles_qR2.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="rightVolumetricBoundaries.nii", intensity_image="qR2_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii right_profiles_qR2.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="leftVolumetricBoundaries.nii", intensity_image="Bieloschowsky_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii left_profiles_Biel.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="rightVolumetricBoundaries.nii", intensity_image="Bieloschowsky_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii right_profiles_Biel.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="leftVolumetricBoundaries.nii", intensity_image="PD_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii left_profiles_PD.nii.gz
rm -R profilesDir/

python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="rightVolumetricBoundaries.nii", intensity_image="PD_resampled.nii.gz", save_data=True, output_dir="profilesDir", file_name="profilesOut.nii")'
mv profilesDir/profilesOut_lps-data.nii right_profiles_PD.nii.gz
rm -R profilesDir/

rm -R rightVolumetric/
rm -R leftVolumetric/
rm -R cruiseRight_brainMask_folder/
rm -R cruiseLeft_brainMask_folder/
rm -R cruiseRight_whiteMatter_folder/
rm -R cruiseLeft_whiteMatter_folder/




