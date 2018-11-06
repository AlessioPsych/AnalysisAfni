#!/bin/bash

SURFACE_LEVELSET_FILE=$1
INTENSITY_FILE=$2

if [ -z "$1" ]
then
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'SURFACE_LEVELSET_FILE=$1, layering ourput'
echo 'INTENSITY_FILE=$2, volume from which extract the profiles' 
echo
exit 1
fi

cp $SURFACE_LEVELSET_FILE _ttt_SURFACE.nii.gz
cp $INTENSITY_FILE _ttt_INTENSIY.nii.gz

run_nighres_command.sh 'laminar.profile_sampling(profile_surface_image="_ttt_SURFACE.nii.gz", intensity_image="_ttt_INTENSIY.nii.gz", save_data=True, output_dir="profiles_output", file_name="profiles")'

rm _ttt_SURFACE.nii.gz
rm _ttt_INTENSIY.nii.gz
#rm _levelset_WM_levelset.nii.gz
#rm _levelset_GM_levelset.nii.gz


# try the topology, seems to be working
# try the sigmoid gradient
# try the resolution


#nighres.laminar.profile_sampling(profile_surface_image, intensity_image, save_data=False, output_dir=None, file_name=None)
