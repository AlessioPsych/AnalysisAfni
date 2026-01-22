#!/bin/bash

docker run -v /media/alessiofracasso/DATADRIVE1/Flanker/derivatives/Freesurfer_output/sub-01/Freesurfer_result/SUMA:/testNighres --rm -it nighres

cd /testNighres

# white matter levelset
python3 -c 'import nighres; surface.probability_to_levelset(probability_image="white_matter_mask.nii.gz", save_data=True, output_dir="whiteMatter", file_name="whiteMatterLevelset.nii.gz")'

cp whiteMatter/*.nii.gz whiteMatterLevelset.nii.gz
