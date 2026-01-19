#!/bin/tcsh -xef

set subj  = $argv[1]
#set subj  = 'sub-01'
set maindir = $PWD
set targetdir = "derivatives/images_processing_afni"
set datadir = "derivatives/processing_afni"

echo $subj
echo $maindir
echo $targetdir
echo $datadir

cd $maindir
cd $datadir
cd {$subj}.results/

@chauffeur_afni                     \
        -ulay    anat_final.{$subj}+tlrc         \
        -prefix  ${maindir}/${targetdir}/{$subj}-anatWarped         \
        -montx 5 -monty 3               \
        -set_xhairs OFF                 \
        -label_mode 1 -label_size 4     \
        -do_clean  


3dTcat -prefix delMeTemp.nii.gz pb02.{$subj}.r01.volreg+tlrc'[1]'

@chauffeur_afni                     \
        -ulay    delMeTemp.nii.gz         \
        -prefix  ${maindir}/${targetdir}/{$subj}-volregWarped         \
        -montx 5 -monty 3               \
        -set_xhairs OFF                 \
        -label_mode 1 -label_size 4     \
        -do_clean  

rm delMeTemp.nii.gz


