#!/bin/bash

set -euo pipefail

# Instructions to run:
# sudo docker run -v /media/alessiofracasso/DATADRIVE1/Flanker:/nighresDataFolder --rm -it nighres
# cd /nighresDataFolder
# bash CD_estimation_volumetricLayering.sh

maindir="/media/alessiofracasso/DATADRIVE1/Flanker/"
targetdir="derivatives/Freesurfer_output"

echo "$maindir"
echo "$targetdir"

# Create subjList.txt file if it does not exist
if [ ! -f subjList.txt ]; then
    ls | grep '^sub-' > subjList.txt
fi

# Loop across participants
while read -r subj; do
    cd "$maindir/$targetdir/$subj/Freesurfer_result/SUMA"

    echo "$PWD"
    echo "$subj"

    # white matter levelset
    python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="white_matter_mask.nii.gz", save_data=True, output_dir="whiteMatter", file_name="whiteMatterLevelset.nii.gz")'

    cp whiteMatter/*.nii.gz whiteMatterLevelset.nii.gz
    chmod 777 whiteMatterLevelset.nii.gz

    # grey matter levelset
    python3 -c 'import nighres; nighres.surface.probability_to_levelset(probability_image="gray_matter_mask_out.nii.gz", save_data=True, output_dir="grayMatter", file_name="grayMatterLevelset.nii.gz")'

    cp grayMatter/*.nii.gz grayMatterLevelset.nii.gz
    chmod 777 grayMatterLevelset.nii.gz

    # volumetric layering
    python3 -c 'import nighres; nighres.laminar.volumetric_layering(inner_levelset="whiteMatterLevelset.nii.gz", outer_levelset="grayMatterLevelset.nii.gz", n_layers=3, topology_lut_dir=None, save_data=True, output_dir="volumetric", file_name="volumetricData.nii.gz")'

    cp volumetric/volumetricData_layering-boundaries.nii.gz volumetricData_layering_boundaries.nii.gz
    cp volumetric/volumetricData_layering-depth.nii.gz volumetricData_layering_depth.nii.gz
    cp volumetric/volumetricData_layering-layers.nii.gz volumetricData_layering_layers.nii.gz
    chmod 777 volumetricData_layering_boundaries.nii.gz
    chmod 777 volumetricData_layering_depth.nii.gz
    chmod 777 volumetricData_layering_layers.nii.gz

    # profile sampling
    python3 -c 'import nighres; nighres.laminar.profile_sampling(profile_surface_image="volumetricData_layering_boundaries.nii.gz", intensity_image="Freesurfer_result_SurfVol.nii", save_data=True, output_dir="profilesT1", file_name="profilesT1.nii.gz")'

    cp profilesT1/*.nii.gz profilesT1.nii.gz
    chmod 777 profilesT1.nii.gz

    # clean up
    rm -r whiteMatter/ grayMatter/ volumetric/ profilesT1/

done < subjList.txt

