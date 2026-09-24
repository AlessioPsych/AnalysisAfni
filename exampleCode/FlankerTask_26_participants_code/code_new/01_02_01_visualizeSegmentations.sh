#!/bin/tcsh -xef

set subj  = $argv[1]
#set subj  = 'sub-01'
set maindir = $PWD
set targetdir = "derivatives/images_processing_Freesurfer"
set datadir = "derivatives/Freesurfer_output"

echo $subj
echo $maindir
echo $targetdir
echo $datadir

cd $maindir
cd $datadir
cd {$subj}
cd Freesurfer_result
cd SUMA

@chauffeur_afni                     \
        -ulay    brainmask.nii.gz \
        -olay aparc+aseg_REN_gm.nii.gz \
        -montx 5 -monty 5               \
        -set_xhairs OFF                 \
        -label_mode 0 -label_size 0     \
        -cbar ROI_i256	\
        -prefix delMe \
        -box_focus_slices brainmask.nii.gz \
        -prefix  ${maindir}/${targetdir}/{$subj}-segmentation         \
        -do_clean  



