#!/bin/bash

WHITE_MATTER_FILE_RAI=$1
CSF_MATTER_PROB_RAI=$2
GRAY_MATTER_PROB_RAI=$3
WHITE_MATTER_PROB_RAI=$4
ANATOMY_RAI=$5

if [ -z "$1" ]
then
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'WHITE_MATTER_FILE_RAI=$1, white matter mask file'
echo 'CSF_MATTER_PROB_RAI=$2'
echo 'GRAY_MATTER_PROB_RAI=$3'
echo 'WHITE_MATTER_PROB_RAI=$4'
echo 'ANATOMY_RAI=$5'
echo
exit 1
fi

cp $WHITE_MATTER_FILE_RAI _ttt_WM.nii.gz
cp $CSF_MATTER_PROB_RAI _ttt_CSF_prob.nii.gz
cp $GRAY_MATTER_PROB_RAI _ttt_GM_prob.nii.gz
cp $WHITE_MATTER_PROB_RAI _ttt_WM_prob.nii.gz
cp $ANATOMY_RAI _ttt_anatomy.nii.gz

3dcalc -a _ttt_CSF_prob.nii.gz -b _ttt_anatomy.nii.gz -expr 'step(not(b))+a' -prefix _ttt_CSF_prob_full.nii.gz

#3dmerge -1blur_fwhm 0.5 -prefix _ttt_CSF_prob_full_blur.nii.gz _ttt_CSF_prob_full.nii.gz
#3dmerge -1blur_fwhm 0.5 -prefix _ttt_GM_prob_blur.nii.gz _ttt_GM_prob.nii.gz
#3dmerge -1blur_fwhm 0.5 -prefix _ttt_WM_prob_blur.nii.gz _ttt_WM_prob.nii.gz

#probability_to_levelset.sh _ttt_WM.nii.gz True $PWD _levelset_WM.nii.gz

#probability_to_levelset.sh _ttt_GM.nii.gz True $PWD _levelset_GM.nii.gz

#3dcalc -a _ttt_WM.nii.gz -b _ttt_GM.nii.gz -c _ttt_ANAT_ss.nii.gz -expr 'and(iszero(step(a+b)),step(c))+iszero(c)' -prefix _ttt_CSF.nii.gz

#3dcalc -a _ttt_WM.nii.gz -b _ttt_GM.nii.gz -expr 'and(iszero(step(a+b)),step(a+b))+iszero(step(a+b))' -prefix _ttt_CSF.nii.gz

#probability_to_levelset.sh _ttt_CSF.nii.gz True $PWD _levelset_CSF.nii.gz

#3dresample -dxyz 0.8 0.8 0.8 -inset _levelset_WM_levelset.nii.gz -rmode Lin -prefix _levelset_WM_levelset_res.nii.gz
#3dresample -dxyz 0.8 0.8 0.8 -inset _levelset_GM_levelset.nii.gz -rmode Lin -prefix _levelset_GM_levelset_res.nii.gz
#3dresample -dxyz 0.8 0.8 0.8 -inset _levelset_CSF_levelset.nii.gz -rmode Lin -prefix _levelset_CSF_levelset_res.nii.gz
#3dresample -master _levelset_WM_levelset_res.nii.gz -inset _ttt_WM.nii.gz -rmode NN -prefix _ttt_WM_res.nii.gz

#levelset_to_probability.sh _levelset_WM_levelset.nii.gz _ttt_WM_PROB.nii.gz $SMOOTH_PAR

#levelset_to_probability.sh _levelset_GM_levelset.nii.gz _ttt_GM_PROB.nii.gz $SMOOTH_PAR

#levelset_to_probability.sh _levelset_CSF_levelset.nii.gz _ttt_CSF_PROB.nii.gz $SMOOTH_PAR
 
run_nighres_command.sh 'cortex.cruise_cortex_extraction(init_image="_ttt_WM.nii.gz", wm_image="_ttt_WM_prob.nii.gz", gm_image="_ttt_GM_prob.nii.gz", csf_image="_ttt_CSF_prob_full.nii.gz",topology="yes",save_data=True, output_dir="cruise_output",file_name="cruise_segmentation")'

rm _ttt_WM.nii.gz
rm _ttt_CSF_prob.nii.gz
rm _ttt_GM_prob.nii.gz
rm _ttt_WM_prob.nii.gz
rm _ttt_*

# try the topology, seems to be working
# try the sigmoid gradient
# try the resolution


