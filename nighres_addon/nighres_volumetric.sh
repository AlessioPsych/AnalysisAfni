#!/bin/bash

WHITE_MATTER_FILE_RAI=$1
GRAY_MATTER_FILE_RAI=$2

if [ -z "$1" ]
then
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'WHITE_MATTER_FILE_RAI=$1, white matter mask file'
echo 'GRAY_MATTER_FILE_RAI=$2, gray matter (ribbon) file' 
echo
exit 1
fi

cp $WHITE_MATTER_FILE_RAI _ttt_WM.nii.gz
cp $GRAY_MATTER_FILE_RAI _ttt_GM.nii.gz

probability_to_levelset.sh _ttt_WM.nii.gz True $PWD _levelset_WM.nii.gz

probability_to_levelset.sh _ttt_GM.nii.gz True $PWD _levelset_GM.nii.gz

run_nighres_command.sh 'laminar.volumetric_layering(inner_levelset="_levelset_WM_levelset.nii.gz", outer_levelset="_levelset_GM_levelset.nii.gz", n_layers=8, save_data=True, output_dir="volumetric_output", file_name="volumetric")'

rm _ttt_WM.nii.gz
rm _ttt_GM.nii.gz
rm _levelset_WM_levelset.nii.gz
rm _levelset_GM_levelset.nii.gz


# try the topology, seems to be working
# try the sigmoid gradient
# try the resolution


